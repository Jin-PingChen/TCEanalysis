# Provide persistent string->string mapping (backed by sqlite file containing a single table). Multiple mappings can be created by passing different identifiers to the constructor.
from builtins import range
from builtins import object
import sys
from sqlalchemy import create_engine, Table, Column, MetaData, sql, String
from sqlalchemy.dialects.sqlite import INTEGER, FLOAT, VARCHAR
from sqlalchemy.orm import declarative_base
from sqlalchemy.orm import sessionmaker
from sqlalchemy.exc import IntegrityError

_debug = False

def get_db(filename):
    db_url = f'sqlite:///{filename}.db'
    db = create_engine(db_url, echo=_debug)
    connection = db.connect()
    return db, connection

md = MetaData()
strings_map = Table("strings_map", md,
                    Column("key", String, primary_key=True),
                    Column("value", String))

Base = declarative_base()
class StringsMap(Base):
    __tablename__ = "strings_map"
    key = Column(String, primary_key=True)
    value = Column(String)


class LocalStringsCache(object):
    def __init__(self, storeName):
        self._db, self._connection = get_db(storeName)
        self._storeName = storeName
        md.create_all(self._db)
        self._session = sessionmaker(bind=self._db)()

    def get_value(self, key):
        results = self._connection.execute(
            sql.select(*(strings_map.c.value,)).select_from(strings_map).where(
                sql.and_(
                    strings_map.c.key == key
                )
            )
        )
        ret = results.fetchall()
        if len(ret) < 1:
            return None
        return ret[0][0]

    def all_matching_values_source(self, prefix):
        results = self._connection.execute(
            sql.select(*(strings_map.c.key, strings_map.c.value,)).select_from(strings_map).where(
                strings_map.c.key.startswith(prefix)
            )
        )
        for u, v in results:
            yield (u, v)

    def insert_value(self, key, value):
        existing_row = self._session.query(StringsMap).filter_by(key=key).first()
        if existing_row:
            existing_row.value = value
        else:
            sm = StringsMap(key=key, value=value)
            self._session.add(sm)
        self._session.commit()

    def update_value(self, key, newValue):
        row = self._session.query(StringsMap).filter_by(key=key).first()
        if row is None:
            raise Exception("Couldn't find StringsMap row matching key={}".format(key))
        row.value = newValue
        self._session.commit()

    def clear_all(self):
        self._session.query(StringsMap).delete()
        self._session.commit()


def testAll():
    sc1 = LocalStringsCache("test_sc1")
    sc1.clear_all()

    sc1.insert_value("test1", "test1")
    sc1.update_value("test1", "----1----")
    sc1.update_value("test1", "----2----")
    sc1.update_value("test1", "----3----")
    sc1.update_value("test1", "test1")

    for i in range(10):
        sc1.insert_value("x.{}".format(i), f"test_seq_{i}")
    for i in range(10):
        sc1.update_value("x.{}".format(i), f"updated_seq_{i}")

    for u, v in sc1.all_matching_values_source("x."):
        print(u, v[:10])

    print(sc1.get_value("test1"))
    print(sc1.get_value("test1X"))
    assert sc1.get_value("test1") == "test1"
    return 0


if __main__ == "__main__":
    import sys
    sys.exit(testAll())
