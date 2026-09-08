# mysql_rnafold.py
from builtins import str
from builtins import object
import sys
from sqlalchemy import create_engine, Table, Column, Integer, Text, String, SmallInteger, Float, MetaData
from sqlalchemy.dialects.mysql import BLOB
from sqlalchemy.orm import declarative_base
from sqlalchemy.orm import sessionmaker
import config


def testHostAvailable(host, port):
    return True


def getMysqlHostFromConnectionString(connString):
    if connString.startswith("sqlite:///"):
        return ("localhost", 0)
    from re import match
    result = match("mysql([^:]+)?://\w+:\w+@([a-zA-z0-9-._]+)(:[0-9]+)?/\w+", connString)
    if not result:
        raise ValueError(f"Invalid connection string: {connString}")
    host = result.group(2)
    portMatch = result.group(3)

    if not portMatch:
        port = 3306
    else:
        assert portMatch[0] == ":"
        port = int(portMatch[1:])

    return (host, port)


_db_cache = {}


def db(config_instance=None):
    if config_instance is None:
        config_instance = config.config.termfold

    cache_key = config_instance.make_sqlite_host_connection() if config_instance.use_sqlite else config_instance.mysql_host_connection
    if cache_key in _db_cache:
        return _db_cache[cache_key]

    if config_instance.use_sqlite:
        conn_str = config_instance.make_sqlite_host_connection()
        _db = create_engine(
            conn_str,
            echo=False,
            pool_recycle=3600,
            connect_args={"check_same_thread": False}
        )
    else:
        if config.run_without_mysql_server:
            return None
        conn_str = config_instance.mysql_host_connection
        _db = create_engine(
            conn_str,
            echo=False,
            echo_pool=True,
            pool_recycle=120
        )

    connection = None
    connectionInfo = getMysqlHostFromConnectionString(conn_str)
    if not config.run_without_mysql_server or config_instance.use_sqlite:
        try:
            connection = _db.connect()
        except Exception as e:
            print("Failed to connect to database on %s" % conn_str)
            print(e)
            sys.exit(-1)

    _db_cache[cache_key] = _db
    return _db


md = MetaData()

sequences2 = Table("sequences2", md,
                  Column("id", Integer, primary_key=True),
                  Column("alphabet", SmallInteger),
                  Column("sequence", BLOB),
                  Column("source", Integer))

sequence_series2 = Table("sequence_series2", md,
                        Column("sequence_id", Integer, primary_key=True),
                        Column("content", BLOB),
                        Column("source", Integer, primary_key=True),
                        Column("ext_index", Integer))

sequence_series2_updates = Table("sequence_series2_updates", md,
                                 Column("dummy_id", Integer, primary_key=True),
                                 Column("sequence_id", Integer),
                                 Column("content", BLOB),
                                 Column("source", Integer),
                                 Column("ext_index", Integer))

sequence_floats2 = Table("sequence_floats2", md,
                         Column("sequence_id", Integer, primary_key=True),
                         Column("value", Float),
                         Column("source", Integer, primary_key=True))


Base = declarative_base()


class Sequence2(Base):
    __tablename__ = "sequences2"
    id = Column(Integer, primary_key=True)
    alphabet = Column(SmallInteger)
    sequence = Column(BLOB)
    source = Column(Integer)


class SequenceSeries2(Base):
    __tablename__ = "sequence_series2"
    sequence_id = Column(Integer, primary_key=True)
    content = Column(BLOB)
    source = Column(Integer, primary_key=True)
    ext_index = Column(Integer)


class SequenceSeries2Updates(Base):
    __tablename__ = "sequence_series2_updates"
    dummy_id = Column(Integer, primary_key=True)
    sequence_id = Column(Integer)
    content = Column(BLOB)
    source = Column(Integer)
    ext_index = Column(Integer)


class SequenceFloats2(Base):
    __tablename__ = "sequence_floats2"
    sequence_id = Column(Integer, primary_key=True)
    value = Column(Float)
    source = Column(Integer, primary_key=True)


engine = db()
if engine:
    Session = sessionmaker(bind=engine)
    Base.metadata.create_all(engine)
    print(f"Successfully connected to database: {config.config.termfold.make_sqlite_host_connection()}")
else:
    raise Exception("Database engine initialization failed! Please check config.py configuration")


class Alphabets(object):
    DNA = 1
    RNA = 2
    RNA_Huff = 3


class Sources(object):
    External = 1
    Computed = 2
    CDSwith3primeFlankingRegion = 3
    CDSwith3primeFlankingRegion_DontExcludeNextORF = 4
    ShuffleCDSv2_matlab = 10
    ShuffleCDSv2_python = 11
    ShuffleCDS_vertical_permutation_1nt = 12
    ShuffleCDS_synon_perm_and_3UTR_nucleotide_permutation = 20
    ShuffleCDS_synon_perm_and_3UTR_nucleotide_permutation_Including_Next_CDS = 21
    ShuffleCDS_synon_perm_and_3UTR_nucleotide_permutation_Including_Next_CDS_Constant_Overlaps = 22
    RNAfoldEnergy_SlidingWindow40_v2 = 102
    RNAfoldEnergy_SlidingWindow40_v2_alt = 103
    RNAfoldEnergy_SlidingWindow40_v2_native_temp = 110
    RNAfoldEnergy_SlidingWindow30_v2 = 113
    RNAfoldEnergy_SlidingWindow50_v2 = 115
    CDS_length_nt = 201
    PA_paxdb_single_assay_or_weighted_average = 202
    GC_content_all_CDS = 203
    GC_content_codon_pos_1 = 204
    GC_content_codon_pos_2 = 205
    GC_content_codon_pos_3 = 206
    MFE_mean_window_40nt_estimated = 207
    GC_content_SlidingWindow40 = 208
    Purine_content_SlidingWindow40 = 209
    StopCodon_content_SlidingWindow40 = 210
    StopCodon_content_SlidingWindow30 = 211
    StopCodon_content_SlidingWindow50 = 212
    TEST_StepFunction_BeginReferenced = 801
    TEST_StepFunction_EndReferenced = 802


windowWidthsByComputationTag = {
    Sources.RNAfoldEnergy_SlidingWindow40_v2: 40,
    Sources.RNAfoldEnergy_SlidingWindow40_v2_native_temp: 40,
    Sources.TEST_StepFunction_BeginReferenced: 40,
    Sources.TEST_StepFunction_EndReferenced: 40,
    Sources.GC_content_SlidingWindow40: 40,
    Sources.Purine_content_SlidingWindow40: 40,
    Sources.StopCodon_content_SlidingWindow40: 40,
    Sources.RNAfoldEnergy_SlidingWindow30_v2: 30,
    Sources.StopCodon_content_SlidingWindow30: 30,
    Sources.RNAfoldEnergy_SlidingWindow50_v2: 50,
    Sources.StopCodon_content_SlidingWindow50: 50
}


def getWindowWidthForComputationTag(computationTag: int) -> int:
    width = windowWidthsByComputationTag.get(computationTag, None)
    if width is None:
        raise ValueError("Unsupported computation-tag {}".format(computationTag))
    else:
        return width
