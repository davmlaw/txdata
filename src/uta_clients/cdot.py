from cdot.hgvs.dataproviders import HGVSDataNotAvailableError, RESTDataProvider


def connect():
    return RESTDataProvider()
