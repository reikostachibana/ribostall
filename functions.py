import ribopy
from ribopy import Ribo
from Fasta import FastaFile
from ribopy.core.get_gadgets import get_region_boundaries, get_reference_names
import pandas as pd

def get_sequence(ribo_object, reference_file_path, alias):
    """
    Retrieves the sequences of transcripts from a reference FASTA file.

    Parameters:
        ribo_object (Ribo): The Ribo object containing ribosome profiling data.
        reference_file_path (str): The file path to the reference FASTA file.
        alias (bool): Whether or not alias is used.

    Returns:
        dict: A dictionary mapping transcript identifiers to their respective sequences.
    """
    transcript_np = ribo_object.transcript_names
    fasta = FastaFile(reference_file_path)
    
    fasta_dict = {e.header: e.sequence for e in fasta}
    if alias == True:
        sequence_dict = {
            transcript.split("|")[4]: fasta_dict[transcript] for transcript in transcript_np
        }
    else:
        sequence_dict = {
        transcript: fasta_dict[transcript] for transcript in transcript_np
    }
    return sequence_dict

def get_psite_offset(ribo_object, exp, mmin, mmax):
    """
    Calculates the P-site offsets for ribosome profiling experiments.

    Parameters:
        ribo_object (Ribo): The Ribo object containing ribosome profiling data.
        exp (str): The name of the experiment for which P-site offsets are calculated.
        mmin (int): The minimum read length considered for P-site offset calculation.
        mmax (int): The maximum read length considered for P-site offset calculation.

    Returns:
        dict: A dictionary mapping transcript identifiers to their respective P-site offsets.
    """
    df = (ribo_object.get_metagene("start", experiments=exp,\
                                   range_lower=mmin, range_upper=mmax,\
                                   sum_lengths=False,\
                                   sum_references=True))

    p_site = {}

    for index, row in df.iterrows():
        max_value_index = row.iloc[35:41].idxmax()
        offset = -1 * max_value_index + 1

        p_site[index[1]] = offset 
    return p_site

def get_cds_range_lookup(ribo_object):
    """
    Create a dict of gene to CDS ranges.
    
    Parameters:
        ribo_object (Ribo): The Ribo object containing ribosome profiling data.

    Returns:
        dict: A dictionary mapping transcript identifiers to the start and stop positions of CDS.
    """
    names = get_reference_names(ribo_object._handle)
    if ribo_object.alias is not None:
        names = map(ribo_object.alias.get_alias, names)
    
    boundaries = get_region_boundaries(ribo_object._handle)
    cds_ranges = [boundary[1] for boundary in boundaries]
    boundary_lookup = dict(zip(list(names), cds_ranges))

    return boundary_lookup