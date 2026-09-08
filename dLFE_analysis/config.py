#!/usr/bin/env python
# coding: utf-8
import logging
from dataclasses import dataclass
import subprocess
import os

# ===================== Utility helper functions (USED in TCE‑Analysis, DO NOT DELETE) =====================
def get_db_connection_string():
    if ConfigTermfold.use_sqlite:
        return ConfigTermfold.make_sqlite_host_connection()
    else:
        return ConfigTermfold.mysql_host_connection


def process_id() -> str:
    hostname = subprocess.check_output("hostname", shell=True).decode().strip()
    pid = os.getpid()
    return f"{hostname}_{pid}"


def init_logging(log_file_name: str = None):
    if log_file_name is None:
        log_file_name = f"process_{process_id()}.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(log_file_name, encoding="utf-8"),
            logging.StreamHandler()
        ]
    )


def ensure_dir_exists(dir_path: str):
    os.makedirs(dir_path, exist_ok=True)

# ===================== Legacy Termfold config: REQUIRED for codon_randomization / mysql_rnafold import compatibility =====================
@dataclass(frozen=True)
class ConfigTermfold:
    use_sqlite = True
    host = "172.**.*.**"
    port = ****
    db = 2
    password = "******"

    sqlite_base_path = "/data/yourname/RTSanalysis/data"
    sqlite_default_db = "rnafold_data"

    @classmethod
    def make_sqlite_host_connection(cls, filename=None):
        if filename is None:
            filename = cls.sqlite_default_db
        return f'sqlite:///{cls.sqlite_base_path}/{filename}.db'

    mysql_host_connection = '***'


# ===================== TCE‑Analysis core config: UTNI / TeRe / Leading / Leading_RTS =====================
@dataclass(frozen=True)
class ConfigTCE:
    BASE_DATA_DIR: str = "/data/chenjp/RTSanalysis/dLFE_code"
    INPUT_GENOME_FOLDER: str = "/data/chenjp/RTSanalysis/dLFE_code/test_genome"

    OUT_FNA_TeRe: str = "/data/chenjp/RTSanalysis/dLFE_code/allfna/TeRe_fna"
    OUT_FNA_LEADING_RTS: str = "/data/chenjp/RTSanalysis/dLFE_code/allfna/Leading_RTS_fna"
    OUT_FNA_UTNI_LEADING: str = "/data/chenjp/RTSanalysis/dLFE_code/allfna/UTNI_Leading_fna"

    TCE_LOG_BASE_NAME: str = "tce_analysis"

    # ---------------------- Gene identification thresholds ----------------------
    UTNI_INTERGENIC_MIN: int = 0
    UTNI_INTERGENIC_MAX: int = 20
    LEADING_LONG_GAP: int = 100   # align with original core_leadinggenes logic

    # ---------------------- ΔLFE calculation parameters ----------------------
    WINDOW_WIDTH: int = 30
    UTR_SPAN: int = 90
    CDS_SPAN: int = 90
    RANDOMIZATION_DEPTH: int = 10
    SEQUENCE_SAMPLING_FRACTION: int = 15

    # ---------------------- Heatmap plotting config ----------------------
    HEATMAP_CMAP_COLORS = ['#2166ac', '#ffff00', '#b10026']
    HEATMAP_N_BINS: int = 50
    HEATMAP_VMIN: float = -2.0
    HEATMAP_VMAX: float = 2.0


# ===================== Top‑level aggregate config object =====================
@dataclass(frozen=True)
class Config:
    termfold = ConfigTermfold()
    tce = ConfigTCE()
    default = termfold


config = Config()

run_without_mysql_server = True
computation_monitor_app = '---pushover‑token---'
computation_monitor_group = '---pushover‑group---'

MatlabPath = ""
codonwBasePath = '/data/chenjp/miniconda3/envs/tceanalysis/bin/codonW'
ENCprimeBasePath = '/data/chenjp/ENCprime‑master/bin/'
