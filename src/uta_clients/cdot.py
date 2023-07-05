from cdot.hgvs.dataproviders import (
    HGVSDataNotAvailableError,
    NotImplementedError,
    RESTDataProvider,
)


def connect():
    return RESTDataProvider()
