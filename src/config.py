from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent

DATA_DIR = PROJECT_ROOT / "data"
KERNEL_DIR = DATA_DIR / "kernels"

STATIONS_CSV = DATA_DIR / "stations.csv" 
PLANTERY_METAKERNEL_TXT = DATA_DIR / "planetaryMetaK.txt"
EARTH_LATEST_HIGH_PRESCISION = KERNEL_DIR / "earth_latest_high_prec.bpc"
EXAMPLE_SPICEKERNEL = KERNEL_DIR / "example.spicekernel"