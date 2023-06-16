import src.uta_clients.json_data_provider as cdot
import hgvs.dataproviders.uta as uta

import datetime
import pytest
from typing import List, Dict

def just_values(dictionaries: List[Dict]) -> List[List]:
    """Like the opposite of dict(), when dict() can't be used.
    (i.e. to compare just the values of a dictionary against a list of values)"""
    return [list(dictionary.values()) for dictionary in dictionaries]

def equal_regardless_of_order(list, other):
    """To compare lists that are neither hashable or stortable."""
    for item in list:
        try:
            other.remove(item)
        except ValueError:
            return False
    return True

def contains(big : list, smol : list) -> bool:
    """Used when either UTA or CDOT has a larger list, and the larger
    list contains everything found in the smaller list."""
    for item in smol:
        if not item in big:
            return False
    return True

#########################################################################
# u = UTA
# c = CDot

@pytest.mark.skip(reason="super slow")
def test_cdot_seq_e():
    """Existing seq."""
    u = uta.connect().get_seq("NC_000007.13")
    c = cdot.connect().get_seq("NC_000007.13")
    assert u == c

@pytest.mark.vcr  
def test_cdot_seq_ne():
    """Nonexisting seq."""
    with pytest.raises(uta.HGVSDataNotAvailableError):
        cdot.connect().get_seq("fake")
    with pytest.raises(uta.HGVSDataNotAvailableError):
        uta.connect().get_seq("fake")

#@pytest.mark.skip(reason="just not fast")
@pytest.mark.vcr
def test_cdot_seq_e_indices():
    """Existing seq with start and end indices."""
    u = uta.connect().get_seq("NC_000007.13", 10000, 10050)
    c = cdot.connect().get_seq("NC_000007.13", 10000, 10050)
    assert u == c

@pytest.mark.skip(reason="cdot has different transcripts than uta")
#@pytest.mark.vcr
def test_cdot_acs_for_protein_seq_e():
    """Existing seq."""
    u = uta.connect().get_acs_for_protein_seq("MRAKWRKKRMRRLKRKRRKMRQRSK")
    c = cdot.connect().get_acs_for_protein_seq("MRAKWRKKRMRRLKRKRRKMRQRSK")
    assert u == c

@pytest.mark.skip(reason="cdot has different transcripts than uta")
#@pytest.mark.vcr    
def test_cdot_acs_for_protein_seq_ne():
    """Nonexisting seq."""
    u = uta.connect().get_acs_for_protein_seq("fake")
    c = cdot.connect().get_acs_for_protein_seq("fake")
    assert u == c
    
@pytest.mark.skip(reason="not implemented in cdot")
#@pytest.mark.vcr
def test_cdot_acs_for_protein_seq_ne_nonalphabetical():
    """Non-alphabetic character seq."""
    with pytest.raises(RuntimeError):
        cdot.connect().get_acs_for_protein_seq("123")
    with pytest.raises(RuntimeError):
        uta.connect().get_acs_for_protein_seq("123")
        
@pytest.mark.skip(reason="cdot has no datetime")
#@pytest.mark.vcr
def test_cdot_gene_info_e():
    """Existing gene."""
    u = dict(uta.connect().get_gene_info("VHL"))
    c = cdot.connect().get_gene_info("VHL")
    del u["summary"]
    del c["summary"]
    c["added"] = datetime.datetime.fromisoformat(c["added"])
    assert u == c
    
@pytest.mark.vcr
def test_cdot_gene_info_ne():
    """Nonexistng gene."""
    u = uta.connect().get_gene_info("VH")
    c = cdot.connect().get_gene_info("VH")
    assert u == c

@pytest.mark.skip(reason="cdot does not include the following fields for transcripts: hgnc, tx_aseq, alt_aseq, test_exon_set_id, aes_exon_set_id, tx_exon_id, alt_exon_id, exon_aln_id")        
#@pytest.mark.vcr
def test_cdot_tx_exons_e():
    """Existing seqs."""
    u = uta.connect().get_tx_exons("NM_199425.2", "NC_000020.10", "splign")
    c = just_values(cdot.connect().get_tx_exons("NM_199425.2", "NC_000020.10", "splign"))
    assert u == c  
    
@pytest.mark.vcr
def test_cdot_tx_exons_ne():
    """Nonexisting seq."""
    with pytest.raises(ValueError):
        cdot.connect().get_tx_exons("NM_199425.2", "fake", "splign")
    with pytest.raises(uta.HGVSDataNotAvailableError):
        uta.connect().get_tx_exons("NM_199425.2", "fake", "splign")
        
@pytest.mark.vcr
def test_cdot_tx_exons_e_params():
    """Existing seqs with missing param"""
    with pytest.raises(TypeError):
        cdot.connect().get_tx_exons("NM_199425.2", "NC_000020.10")
    with pytest.raises(TypeError):
        uta.connect().get_tx_exons("NM_199425.2", "NC_000020.10")
        
@pytest.mark.vcr
def test_cdot_tx_exons_ne_params():
    """Existing seqs with missing param"""
    with pytest.raises(TypeError):
        cdot.connect().get_tx_exons("NM_199425.2", "fake")
    with pytest.raises(TypeError):
        uta.connect().get_tx_exons("NM_199425.2", "fake")

@pytest.mark.skip(reason="cdot and uta return two different lists of transcript info with some overlap.")  
#@pytest.mark.vcr
def test_cdot_tx_for_gene_e():
    """Existing gene."""
    u = uta.connect().get_tx_for_gene("VHL")
    c = cdot.connect().get_tx_for_gene("VHL")
    assert equal_regardless_of_order(u, c)

@pytest.mark.vcr
def test_cdot_tx_for_gene_ne():
    """Nonexisitng gene."""
    u = uta.connect().get_tx_for_gene("VH")
    c = cdot.connect().get_tx_for_gene("VH")
    assert u == c
    
@pytest.mark.vcr
def test_cdot_tx_for_region_e(): # Need example
    """Existing region."""
    u = uta.connect().get_tx_for_region("NC_000007.13", "splign", 0, 50)
    c = cdot.connect().get_tx_for_region("NC_000007.13", "splign", 0, 50)
    assert u == c
    
@pytest.mark.vcr
def test_cdot_tx_for_region_ne():
    """Nonexisting region"""
    u = uta.connect().get_tx_for_region("fake", "splign", 0, 50)
    c = cdot.connect().get_tx_for_region("fake", "splign", 0, 50)
    assert u == c
    
@pytest.mark.vcr
def test_cdot_tx_for_region_e_params():
    """Existing region with only some params."""
    with pytest.raises(TypeError):
        cdot.connect().get_tx_for_region("NC_000007.13", "splign")
    with pytest.raises(TypeError):
        uta.connect().get_tx_for_region("NC_000007.13", "splign")
        
@pytest.mark.vcr
def test_cdot_tx_for_region_ne_params():
    """Existing region with only some params."""
    with pytest.raises(TypeError):
        cdot.connect().get_tx_for_region("fake", "splign")
    with pytest.raises(TypeError):
        uta.connect().get_tx_for_region("fake", "splign")
    
@pytest.mark.skip(reason="not implemented in cdot")
#@pytest.mark.vcr
def test_cdot_alignments_for_region_e(): # Need example
    """Existing region."""
    u = uta.connect().get_alignments_for_region("NC_000007.13", 0, 50)
    c = cdot.connect().get_alignments_for_region("NC_000007.13", 0, 50)
    assert u == c
    
@pytest.mark.skip(reason="not implemented in cdot")
#@pytest.mark.vcr
def test_cdot_alignments_for_region_ne():
    """Nonexisting region"""
    u = uta.connect().get_alignments_for_region("fake", 0, 50)
    c = cdot.connect().get_alignments_for_region("fake", 0, 50)
    assert u == c
    
@pytest.mark.skip(reason="not implemented in cdot")
#@pytest.mark.vcr
def test_cdot_alignments_region_e_params():
    """Existing region with only some params."""
    with pytest.raises(TypeError):
        cdot.connect().get_alignments_for_region("NC_000007.13")
    with pytest.raises(TypeError):
        uta.connect().get_alignments_for_region("NC_000007.13")
        
@pytest.mark.skip(reason="not implemented in cdot")
#@pytest.mark.vcr
def test_cdot_alignments_region_ne_params():
    """Existing region with only some params."""
    with pytest.raises(TypeError):
        cdot.connect().get_alignments_for_region("fake")
    with pytest.raises(TypeError):
        uta.connect().get_alignments_for_region("fake")
    
@pytest.mark.vcr
def test_cdot_tx_identity_info_e():
    """Existing transcript."""
    u = dict(uta.connect().get_tx_identity_info("NM_199425.2"))
    c = cdot.connect().get_tx_identity_info("NM_199425.2")
    assert u == c
    
@pytest.mark.vcr
def test_cdot_tx_identity_info_ne():
    """Nonexisting transcript."""
    assert cdot.connect().get_tx_identity_info("fake") == None
    with pytest.raises(uta.HGVSDataNotAvailableError):
        uta.connect().get_tx_identity_info("fake")

@pytest.mark.vcr
def test_cdot_tx_info_e():
    """Existing seqs."""
    u = dict(uta.connect().get_tx_info("NM_199425.2", "NC_000020.10", "splign"))
    c = cdot.connect().get_tx_info("NM_199425.2", "NC_000020.10", "splign")
    assert u == c
    
@pytest.mark.xfail
@pytest.mark.vcr
def test_cdot_tx_info_ne():
    """Nonexisting seq."""
    assert cdot.connect().get_tx_info("NM_199425.2", "fake", "splign") == None
    with pytest.raises(uta.HGVSDataNotAvailableError):
        uta.connect().get_tx_info("NM_199425.2", "fake", "splign")
        
@pytest.mark.vcr
def test_cdot_tx_info_e_param():
    """Existing seqs."""
    with pytest.raises(TypeError):
        cdot.connect().get_tx_info("NM_199425.2", "NC_000020.10")
    with pytest.raises(TypeError):
        uta.connect().get_tx_info("NM_199425.2", "NC_000020.10")   
    
@pytest.mark.vcr
def test_cdot_tx_info_ne_param():
    """Existing seqs."""  
    with pytest.raises(TypeError):
        cdot.connect().get_tx_info("NM_199425.2", "fake")
    with pytest.raises(TypeError):
        uta.connect().get_tx_info("NM_199425.2", "fake")    
        
@pytest.mark.vcr
def test_cdot_tx_mapping_options_e():
    """Existing transcript."""
    u = uta.connect().get_tx_mapping_options("NM_000051.3")
    c = just_values(cdot.connect().get_tx_mapping_options("NM_000051.3"))
    #assert u == c
    assert contains(u, c)
       
@pytest.mark.vcr
def test_cdot_tx_mapping_options_ne():
    """Nonexisting transcript."""
    u = uta.connect().get_tx_mapping_options("fake")
    c = just_values(cdot.connect().get_tx_mapping_options("fake"))
    assert u == c
    
@pytest.mark.skip("not implemented in cdot")
#@pytest.mark.vcr
def test_cdot_similar_transcripts_e():
    """Existing transcript."""
    u = uta.connect().get_similar_transcripts("NM_000051.3")
    c = just_values(cdot.connect().get_similar_transcripts("NM_000051.3"))
    assert u == c

@pytest.mark.skip("not implemented in cdot")    
#@pytest.mark.vcr
def test_cdot_similar_transcripts_ne():
    """Nonexisting transcript."""
    u = uta.connect().get_similar_transcripts("fake")
    c = just_values(cdot.connect().get_similar_transcripts("fake"))
    assert u == c
    
@pytest.mark.vcr
def test_cdot_pro_ac_for_tx_ac_e():
    """Existing transcript."""
    u = uta.connect().get_pro_ac_for_tx_ac("NM_000051.3")
    c = cdot.connect().get_pro_ac_for_tx_ac("NM_000051.3")
    assert u == c
    
@pytest.mark.vcr
def test_cdot_pro_ac_for_tx_ac_ne():
    """Nonexisting transcript."""   
    u = uta.connect().get_pro_ac_for_tx_ac("fake")
    c = cdot.connect().get_pro_ac_for_tx_ac("fake")
    assert u == c
    
@pytest.mark.vcr
def test_cdot_assembly_map_e():
    """Existing assembly name."""
    u = dict(uta.connect().get_assembly_map("GRCh38"))
    c = cdot.connect().get_assembly_map("GRCh38")
    assert u == c
    
@pytest.mark.vcr
def test_cdot_assembly_map_ne():
    """Nonexisting assembly name."""
    with pytest.raises(Exception):
        cdot.connect().get_assembly_map("GROUCH")
    with pytest.raises(Exception):
        uta.connect().get_assembly_map("GROUCH")

@pytest.mark.skip(reason="different ways of representing versions")     
#@pytest.mark.vcr
def test_cdot_data_version():
    u = uta.connect().data_version()
    c = cdot.connect().data_version()
    assert u == c

@pytest.mark.skip(reason="different ways of representing versions")   
#@pytest.mark.vcr
def test_cdot_schema_version():
    u = uta.connect().schema_version()
    c = cdot.connect().schema_version()
    assert u == c

@pytest.mark.skip(reason="different ways of representing seq source")  
#@pytest.mark.vcr
def test_cdot_sequence_source():
    u = uta.connect().sequence_source()
    c = cdot.connect().sequence_source()
    assert u == c