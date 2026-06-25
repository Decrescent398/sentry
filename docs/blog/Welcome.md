# Rebuilding Sentry-II From Scratch: Orbit Determination 

What the hell is Sentry-II?. It is NASA's planetary defense system watches every known near-Earth asteroid, fits an orbit to the telescope data, samples tens of thousands of virtual clones across the uncertainty region, propagates them forward a century, and tells you the probability that any of them hit Earth.

Naively, the first thing I did was look up a tutorial on youtube. Not a bad idea in retrospect, but not the best use of time either. I went through Dr. Thomas Albin's playlist on [NEOs](https://www.youtube.com/playlist?list=PLNvIBWkEdZ2hL5be8mQdpTU3BjhKIhD6L) and it was a pretty good introduction and gave me a high level understanding of what was going on. At the time, I felt like I was coding the pinnacle of astrodynamics, but realized that I was in fact just copying code while not really knowing what was going, what parameter impacted others, how were they connected and why it mattered. This eventually led to a lot of technical debt and, I abandoned the project for a couple of months.

I wanted to understand what does it actually mean to _know_ where an asteroid is? What confidence do we have that the math we're doing is correct? 

I figured out that the only way I would truly understand what was going on was by re-implementing the [paper](https://iopscience.iop.org/article/10.3847/1538-3881/ac193f). So I did, using real telescope observations of 433 Eros from the Minor Planet Center, validated against JPL Horizons to within 0.5% of their calculations. This post is about the conceptual surprises where my mental model of the problem was wrong in ways that took weeks to untangle. The code is mostly incidental.

## What the Problem Actually Is

The thing that took the longest to fully internalize, was that a telescope gives you nothing but a direction.

Every observation is just Right Ascension $\alpha$ and Declination $\delta$. No range, no range-rate, no depth, the asteroid could be at **0.3 AU** or **3 AU** and the angles might be identical. So how do you reconstruct a six-dimensional state vector with three position coordinates, and three velocity components from a sequence of two-dimensional measurements?

You use physics (Thank Newton).

If the asteroid obeys Newton's laws (it does), then its trajectory is not arbitrary. The sequence of angles it traces across the sky over time is tightly constrained by the equations of motion. Specifically, the trajectory must satisfy:

$$\ddot{\vec{r}} = -\frac{\mu_\odot}{|\vec{r}|^3}\vec{r} + \sum_i \vec{a}_i(\vec{r}, t)$$

where the first term is the Sun's gravity and the second accumulates perturbations from everything else. The constraint that the trajectory obeys this equation continuously is what makes a sequence of angular measurements informative about a three-dimensional orbit.

The orbit determination problem is then about find the initial state $x_0 = [\vec{r}_0, \dot{\vec{r}}_0]$ such that, when propagated forward through the equations of motion and projected onto the sky plane, the predicted angles match the observed angles as well as possible. Formally, this is a nonlinear least-squares problem defined as:

$$\hat{x} = \arg\min_x \sum_k \frac{(\alpha_k^{\text{pred}}(x) - \alpha_k^{\text{obs}})^2}{\sigma_{\alpha,k}^2} + \frac{(\delta_k^{\text{pred}}(x) - \delta_k^{\text{obs}})^2}{\sigma_{\delta,k}^2}$$

But saying "nonlinear least-squares" massively undersells the complexity. The function $\alpha_k^{\text{pred}}(x)$ requires integrating the equations of motion from epoch to observation time $t_k$, accounting for the position of the observer on the rotating Earth at $t_k$, correcting for the fact that the light you are observing left the asteroid 10 minutes ago, and projecting the resulting geometric vector onto sky-plane angles. Each of these is its own sub-problem with its own failure modes.

## The Solar System Is Not Optional

My first instinct was to start with two-body dynamics with just the Sun and add perturbations later. This is how most textbooks sequence the problem. Boy, was this the wrong instinct.

Over a 90-day arc, the difference between a two-body orbit and a fully perturbed orbit for Eros is not a small residual you can ignore at first pass. Jupiter alone induces positional errors of thousands of kilometers. When you are trying to match sky positions to within an arcsecond which at Eros's distance of ~1.5 AU corresponds to about 1,000 km starting with two-body dynamics means your residuals never drop below the noise floor regardless of how well you tune the estimator. You would be fitting the wrong trajectory, and you would not know it because the residuals would plateau at a "seems wrong but not obviously wrong" level.

So, the integrator includes all eight planets, the Moon, and the eighteen most massive main-belt asteroids. The right way to write the perturber contribution is not just "add the gravitational pull of Jupiter on Eros." You also have to subtract the gravitational pull of Jupiter on the Sun, because you are integrating in a heliocentric frame and the Sun is not inertially fixed when Jupiter is nearby. This indirect term is easy to miss and, if missed, introduces a spurious force equal to the Sun's acceleration from the perturbers, a force that does not actually exist from Eros's perspective. The full equation is:

$$\ddot{\vec{r}} = -\frac{\mu_\odot}{|\vec{r}|^3}\vec{r} - \sum_i \mu_i \left(\frac{\vec{r} - \vec{r}_i}{|\vec{r} - \vec{r}_i|^3} + \frac{\vec{r}_i}{|\vec{r}_i|^3}\right)$$

The second term inside the brackets is a correction to the coordinate system, not an additional force on the asteroid. I initially implemented it without the indirect term. The orbit would not converge below 3 arcseconds RMS no matter what I did, and for a week I blamed the differential corrector.

Python

```
def n_body_ode(t, state):
    r, v = state[:3], state[3:]
    a = -MU_SUN * r / np.linalg.norm(r)**3

    for name, mu in PERTURBERS:
        r_body, _ = spice.spkezr(name, t, 'J2000', 'NONE', 'SUN')
        r_body = r_body[:3]
        diff = r - r_body
        a += -mu * (diff / np.linalg.norm(diff)**3       # direct
                  + r_body / np.linalg.norm(r_body)**3)  # indirect

    return np.concatenate([v, a])
```

The DOP853 integrator at tolerance $10^{-12}$ keeps position errors below a few meters over 90 days. The tolerance sounds aggressive until you remember that you are matching against observations sensitive to ~100 km accuracy, and numerical integration error is about the only error in the pipeline you can actually drive to zero.

## SPICE: The Infrastructure That Bit Me

I knew going in that SPICE would require some setup. I did not expect it to consume more debugging time than the orbital mechanics.

SPICE (Spacecraft Planet Instrument Camera-matrix Events) is NASA's toolkit for computing precise positions of solar system bodies, handling coordinate frame transformations, and doing time-scale conversions. It is what JPL actually uses internally, which means any discrepancy between my planetary positions and theirs is immediately interpretable.

The modularity that makes SPICE powerful is also what makes it bite you. Each type of information lives in a separate kernel file: leap seconds (`.tls`), planetary ephemerides (`.bsp`), Earth orientation (`.bpc`), gravitational parameters (`.tpc`), and frame definitions (`.tf`). 

The specific error I hit was loading `codes_300ast_20100725.bsp`, a three-hundred-asteroid ephemeris compiled by Jim Baer. This file uses a custom frame called `ECLIPJ2000_DE405` (frame ID 1900017), which is an ecliptic frame defined relative to the DE405 dynamical ecliptic rather than the modern IAU standard. The frame definition lives nowhere in the standard kernel set. I had to inject it manually into SPICE's kernel pool:

Python

```
spice.load_kernel_data("""
\begindata
FRAME_ECLIPJ2000_DE405       = 1900017
FRAME_1900017_CLASS          = 4
TKFRAME_1900017_SPEC         = 'ANGLES'
TKFRAME_1900017_RELATIVE     = 'J2000'
TKFRAME_1900017_ANGLES       = ( 0.0, 0.0, -84381.412 )
TKFRAME_1900017_AXES         = ( 1, 3, 1 )
TKFRAME_1900017_UNITS        = 'ARCSECONDS'
\begintext
""")
```

The rotation angle, **84381.412 arcseconds**, is the obliquity of the ecliptic **23.44 degrees** converting between the ecliptic and equatorial reference planes. The lesson is that in SPICE, the distinction between "loading data" and "loading the geometry to interpret the data" is entirely your responsibility. The toolkit cannot infer a frame definition from the data that references it.

The more confusing issue, and actually a product of bad implementation was time. SPICE measures time in ephemeris seconds past J2000, also called TDB (Barycentric Dynamical Time) or ET. MPC observation timestamps are ISO 8601 strings in UTC. The MPC orbital epoch is a Modified Julian Date. These are not interchangeable at all. UTC incorporates leap seconds, TDB is a relativistic timescale tied to the solar system barycenter, and MJD and JD differ by a fixed constant (**2400000.5 days**). Getting any conversion wrong by one second produces a position error of ~25 km for Eros, since that is roughly its orbital speed. These errors never raise exceptions, they just produce wrong numbers.

## The Observation Model: Projecting a 6D Trajectory onto the Sky

For each observation, the pipeline needs to predict $(\alpha_k^{\text{pred}}, \delta_k^{\text{pred}})$ given the current state estimate. This projection involves three nested corrections, each of which has a conceptual subtlety.

**Observer position:** The telescope is not at the solar system center, or at Earth's center, or even at Earth's surface center. It is at a specific geodetic latitude, longitude, and altitude. The MPC publishes observatory parallax constants for each terrestrial station two numbers that encode the observer's offset from Earth's rotation axis and equatorial plane. (Space-based observatories don't provide simple geocentric coordinates in the MPC data and require querying entirely separate satellite ephemerides, which is another trap waiting to be sprung). Converting terrestrial constants to an inertial J2000 position at the observation epoch requires rotating through the Earth's orientation at that moment, including precession, nutation, and polar wander. SPICE handles this, but you have to ask it for the right transformation.

**Light-time correction:** The asteroid in the image is not where the asteroid is right now. The light you are looking at left the asteroid $\Delta t = |\vec{r}_\text{ast}(t_\text{emit}) - \vec{r}_\text{obs}(t_\text{obs})| / c$ seconds ago, where $t_\text{obs}$ is the observation epoch. This is roughly 10 minutes for Eros at typical distances. In that time the asteroid has moved ~15,000 km. The correction is iterative: start with $t_\text{emit} = t_\text{obs}$, compute the geometric distance, compute the implied travel time, step back by that amount, get the new position, repeat. Three iterations converge to machine precision.

Python

```
t_emit = t_obs
for _ in range(3):
    r_ast = trajectory(t_emit)[:3]
    delta = r_ast - r_obs
    t_emit = t_obs - np.linalg.norm(delta) / C_KM_S
```

What you actually want is the position of the asteroid at emission, as seen from the observer at reception. Ignoring it gives systematic residuals at the 5-arcsecond level, which completely swamps the observational noise.

**Projection to angles:** Once you have the corrected geometric look-vector $\vec{\rho} = \vec{r}_\text{ast}(t_\text{emit}) - \vec{r}_\text{obs}(t_\text{obs})$, the predicted angles are:

$$\alpha^{\text{pred}} = \text{atan2}(\rho_y, \rho_x), \quad \delta^{\text{pred}} = \arcsin\left(\frac{\rho_z}{|\vec{\rho}|}\right)$$

The residuals in the cost function are $\Delta\alpha\cos\delta$ and $\Delta\delta$. The $\cos\delta$ factor converts RA differences to great-circle distances on the sky, since RA lines of constant value converge at the poles.

**Observation weighting:** Not all telescope observations are equal, an observation from 1998 reduced against USNO-A2.0 carries a typical uncertainty around 0.65 arcseconds. The same station today reduced against Gaia DR3 carries ~0.12 arcseconds. Treating all observations as equal dramatically changes which epochs dominate the orbit solution. I used the Vereš-Farnocchia-Chesley-Chambers (VFCC) table, derived from Chesley, Baer & Monet (2010), which lists per-(station, catalog) $\sigma_\alpha$ and $\sigma_\delta$ from statistical analysis of millions of historical asteroid residuals (debiasing these catalog positions is an entire sub-field of its own). The weight matrix $W$ is diagonal with entries $1/\sigma^2$.

## Differential Correction: What the Math Is Actually Doing

The orbit estimator is Gauss-Newton iteration with Levenberg-Marquardt damping. I want to say something about why it works, not just describe the mechanics.

The idea is: you have a function $f: \mathbb{R}^6 \to \mathbb{R}^{2N}$ mapping an initial state to a vector of RA/Dec predictions across all $N$ observations. You want to find the zero of $e = f(x) - y^{\text{obs}}$. Newton's method would require the Hessian of the cost. Gauss-Newton avoids this by linearizing $f$ around the current estimate and solving a linear system for the step. The linearization is the design matrix $B$:

$$B_{ij} = \frac{\partial e_i}{\partial x_j}\bigg|_{x = x_k}$$

which I compute numerically: nudge each of the six state components by a small $\delta_j$, rerun the full N-body integration, and compute how the residuals change:

$$B_{ij} \approx \frac{e_i(x + \delta_j \hat{e}_j) - e_i(x)}{\delta_j}$$

Each column of $B$ requires one full integration pass. Six columns, plus the baseline, means seven integrations per iteration. The update step minimizes the linearized residual norm:

$$\Delta x = -(B^\top W B + \lambda I)^{-1} B^\top W e$$

The $\lambda I$ Tikhonov term (Levenberg-Marquardt) regularizes the normal equations during early iterations when the linearization is unreliable and the information matrix $C = B^\top W B$ may be nearly singular.

The thing I got wrong conceptually was thinking that the perturbation step size $\delta_j$ is a numerical detail. It is not. It is a physics question.

The six state components live in different physical units and have different sensitivities. A position nudge of 1 km produces a residual change of order $1/|\vec{\rho}| \approx 10^{-9}$ radians per km at 1 AU. A velocity nudge of $10^{-3}$ km/s accumulated over a 90-day arc (roughly $8 \times 10^6$ seconds) becomes a position offset of $8 \times 10^3$ km, which is eight thousand times larger than the position nudge. If you use the same $\delta$ for position and velocity, the velocity columns of $B$ dominate completely. The information matrix looks rank-deficient to the position components. It _is_ effectively rank-deficient, and the normal equations either fail or produce a solution that trusts the velocity observations far more than the position observations, which is wrong.

The fix is to match $\delta_j$ to the physical scale each parameter actually influences in the residuals. I used 1 km for positions and $10^{-4}$ km/s for velocities, chosen so that each nudge produces a similar-sized residual perturbation. This is dimensional analysis, not numerical tuning.

The convergence is fast once the conditioning is right. Four or five iterations reduce the RMS by a factor of 50 or more from a reasonable initial guess. Past that, the iterations are essentially doing nothing; you are in the quadratic convergence regime of Newton's method and the linearization is accurate.

## The Residual RMS: A Single Number That Told Me Everything

I want to make a case for obsessing over the sky-plane RMS residual as a debugging instrument.

The RMS is the square root of $Q / 2N$ the mean squared residual across all observations, in arcseconds. It is a scalar. It is easy to compute. And it has a ground truth: the best ground-based astrometry with Gaia-calibrated star catalogs achieves a residual floor around 0.15 arcseconds from atmospheric turbulence. Any RMS well above this value means something is wrong with the model, and how far above it tells you _how_ wrong.

During development, I kept a log of RMS values after each significant fix:

|**What I did**|**RMS**|
|---|---|
|Initial integration with indirect term missing|148 arcsec|
|Added indirect term|18.4 arcsec|
|Fixed perturbation scale mismatch in $B$|3.7 arcsec|
|Fixed light-time correction|0.91 arcsec|
|Fixed observer geodetic → J2000 rotation|0.31 arcsec|
|Switched from uniform to VFCC weights|0.18 arcsec|
|Converged solution|0.13 arcsec|

Each row is a conceptual error in the physical model, not a code bug in the conventional sense. "Missing indirect term" is a wrong understanding of what heliocentric means when the heliocenter is accelerating. "Wrong perturbation scale" is a wrong understanding of what information the design matrix needs to encode. The residual RMS diagnosed all of them, even when I did not know what to look for.

The median of the converged residuals is 0.07 arcseconds, below the atmospheric floor, which means the observations are a mix of modern Gaia-era measurements (low noise, high weight) and legacy observations (higher noise, lower weight, contributing less to the RMS).

## The Coordinate Frame Problem: When the Right Answer Looks Wrong

When the filter finally converged, I compared my state vector to JPL Horizons. The Y component of my solution was positive. JPL's Y was negative. My Z component was around $+2 \times 10^7$ km. JPL's Z was $-3 \times 10^7$ km. Several velocity components had the wrong sign.

My first instinct was that something fundamental was broken. My second instinct was to check the residuals, they were still 0.13 arcseconds. The orbit was fitting the telescope data beautifully. How do you simultaneously fit 18,000 real observations to sub-arcsecond precision and have the wrong trajectory?

The answer, embarrassingly in retrospect, is that the same physical orbit can be represented by completely different coordinate values in different reference frames. I was comparing my Cartesian state vector (expressed in J2000 equatorial coordinates, with the reference plane being Earth's mean equatorial plane at the J2000 epoch) to the JPL output (which defaulted to ECLIPJ2000 coordinates, with the reference plane being the ecliptic). These two frames are tilted by the obliquity of the ecliptic: **23.439°**.

A rotation of 23.4 degrees about the x-axis maps between them:

$$\vec{r}_\text{equatorial} = R_x(-\varepsilon) \cdot \vec{r}_\text{ecliptic}$$

After applying this rotation, the position components agreed to within **1%** and the velocity components agreed to within **0.5%** the remaining difference being attributable to JPL's multi-decade, multi-mission orbit solution versus my 90-day arc.

## What the Converged Solution Gives You

After convergence, the Gauss-Newton process has produced two things: a best-fit state vector $\hat{x}$ and an information matrix $C = B^\top W B$. The covariance is the inverse:

$$\Sigma = C^{-1} = (B^\top W B)^{-1}$$

This is a $6 \times 6$ matrix encoding the uncertainty in each state component and the correlations between them. The diagonal entries give the formal 1-sigma uncertainties. For Eros with a 90-day arc, the formal position uncertainty is tens of kilometers and the velocity uncertainty is fractions of a mm/s.

This covariance is the output of orbit determination. A state without a covariance is a number without an error bar. What it means to "know" where an asteroid is turns out to mean: knowing the six-dimensional probability distribution over states consistent with the observations. That distribution, under the assumptions of the Gauss-Newton model, is a Gaussian centered at $\hat{x}$ with covariance $\Sigma$.

## Where This Is Going

Orbit determination is the first problem in a longer chain.

The covariance $\Sigma$ defines a six-dimensional ellipsoid of uncertainty in state space, the set of initial conditions that fit the observations within their errors. To compute whether an asteroid might hit Earth, you need to know what that ellipsoid looks like in 100 years. Propagating the covariance forward linearly, $\Sigma(t) = \Phi(t) \Sigma \Phi(t)^\top$, works for short time horizons but breaks down when the orbit passes near a resonance or makes a close approach: the uncertainty distribution becomes non-Gaussian and the linear approximation fails.

The solution, which is what Sentry-II does, is Monte Carlo: sample thousands of virtual asteroids from the uncertainty distribution, propagate each one independently, and let the distribution of outcomes speak for itself.

Python

```
clones = np.random.multivariate_normal(x_hat, Sigma, size=10_000)
```

Ten thousand initial conditions, all statistically consistent with the observed telescope data. Propagating each one forward 100 years through the full N-body integrator is computationally heavy, but fortunately, it is an embarrassingly parallel problem—perfect for mapping horizontally across a cloud-based container infrastructure. Some of them will have close approaches to Earth. A fraction of those will have minimum orbital intersection distances inside Earth's gravitational capture cross-section. That fraction is the impact probability.

The picture is even richer than this, because the uncertainty ellipsoid is not static. As new observations come in, the ellipsoid shrinks. The information matrix accumulates $C \leftarrow C + b_k^\top w_k b_k$ for each new observation and the covariance tightens. An asteroid with a non-negligible impact probability can often be rescued from ambiguity by a few targeted observations that pin down the uncertainty dimension most relevant to the close approach geometry. This is how planetary defense actually works in practice: not just cataloguing, but targeted follow-up observation to resolve the probability distribution before the encounter.

Beyond the Monte Carlo step, there are physical effects I have not yet modeled. The Yarkovsky effect—which is thermal radiation pressure from an asteroid's asymmetric heating and cooling—can shift an orbit by enough to move a near-miss into an impact trajectory over decades. General relativistic corrections add roughly 40 meters per orbit for Eros, which accumulates over centuries. Neither of these is hard to add to the force model once the basic infrastructure is solid.

The long-term goal is an open planetary defense platform. Orbit determination, covariance propagation, virtual asteroid sampling, impact probability estimation—all of it, open-source, documented, reproducible. Sentry-II produces a number. I want a system where you can trace that number backwards through every assumption and modeling choice, inspect the sensitivity to each one, and understand exactly what "**0.0014%** probability of impact in 2078" actually means and does not mean.

Building it is slower than I expected, mostly because each step reveals physical subtleties I thought I understood and did not. The indirect term in the N-body ODE. The frame definitions in SPICE kernels. The timescale conversions. The meaning of a Cartesian state vector without a reference frame. None of these are in the textbooks in the form you need them when something is silently wrong.

But the residual RMS keeps dropping. That is progress.