# uta-clients

A package containing clients of [UTA](https://github.com/biocommons/uta). 

Each acts as a data provider for hgvs, a library of sequence manipulation tools. Ensuring a data provider for hgvs allows users to normalize, validate, and map sequence variants, among other functionalities such as parsing and formatting that are included in the [hgvs](https://github.com/biocommons/hgvs) library.

## Clients

Supported clients include: 
- [uta](https://github.com/biocommons/hgvs/blob/main/src/hgvs/dataproviders/uta.py)
- uta ([rest](https://github.com/ccaitlingo/uta-rest))
- cdot ([rest](https://github.com/SACGF/cdot/blob/main/cdot/hgvs/dataproviders/json_data_provider.py))


## Using with hgvs

    $ pip install hgvs
    $ pip install utaclients

In python:

    >>> import hgvs
    >>> import utaclients
    >>> hdp = ...
    
Assign any of the following supported clients to your hgvs data provider.

    >>> hdp = utaclients.uta.connect()
    >>> hdp = utaclients.uta_rest.connect()
    >>> hdp = utaclients.cdot.connect()

Then you can get starting using all the functionalities of hgvs.

    >>> hdp = utaclients.uta.connect()

    >>> am = hgvs.assemblymapper.AssemblyMapper(hdp, 
    ...     assembly_name='GRCh37', alt_aln_method='splign',
    ...     replace_reference=True)