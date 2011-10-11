"""Perform realignment of BAM files around indels using the GATK toolkit.
"""
import os
import shutil
from contextlib import closing

import pysam

from bcbio import broad
from bcbio.pipeline import log
from bcbio.utils import curdir_tmpdir, file_transaction
from bcbio.distributed.split import parallel_split_combine
from bcbio.pipeline.shared import (split_bam_by_chromosome, configured_ref_file)

# ## Realignment runners with GATK specific arguments

def gatk_realigner_targets(runner, align_bam, ref_file, dbsnp=None,
                           region=None, out_file=None, deep_coverage=False):
    """Generate a list of interval regions for realignment around indels.
    """
    if out_file:
        out_file = "%s.intervals" % os.path.splitext(out_file)[0]
    else:
        out_file = "%s-realign.intervals" % os.path.splitext(align_bam)[0]
    params = ["-T", "RealignerTargetCreator",
              "-I", align_bam,
              "-R", ref_file,
              "-o", out_file,
              "-l", "INFO",
              ]
    if region:
        params += ["-L", region]
    if dbsnp:
        params += ["--known", dbsnp]
    if deep_coverage:
        params += ["--mismatchFraction", "0.30",
                   "--maxIntervalSize", "650"]
    if not (os.path.exists(out_file) and os.path.getsize(out_file) > 0):
        with file_transaction(out_file):
            runner.run_gatk(params)
    return out_file

def gatk_indel_realignment(runner, align_bam, ref_file, intervals,
                           region=None, out_file=None, deep_coverage=False):
    """Perform realignment of BAM file in specified regions
    """
    if out_file is None:
        out_file = "%s-realign.bam" % os.path.splitext(align_bam)[0]
    params = ["-T", "IndelRealigner",
              "-I", align_bam,
              "-R", ref_file,
              "-targetIntervals", intervals,
              "-o", out_file,
              "-l", "INFO",
              ]
    if region:
        params += ["-L", region]
    if deep_coverage:
        params += ["--maxReadsInMemory", "300000",
                   "--maxReadsForRealignment", str(int(5e5)),
                   "--maxReadsForConsensuses", "500",
                   "--maxConsensuses", "100"]
    if not (os.path.exists(out_file) and os.path.getsize(out_file) > 0):
        with curdir_tmpdir() as tmp_dir:
            with file_transaction(out_file):
                runner.run_gatk(params, tmp_dir)
    return out_file

def gatk_realigner(align_bam, ref_file, config, dbsnp=None, region=None,
                   out_file=None, deep_coverage=False):
    """Realign a BAM file around indels using GATK, returning sorted BAM.
    """
    runner = broad.runner_from_config(config)
    runner.run_fn("picard_index", align_bam)
    runner.run_fn("picard_index_ref", ref_file)
    if not os.path.exists("%s.fai" % ref_file):
        pysam.faidx(ref_file)
    if _has_aligned_reads(align_bam, region):
        realign_target_file = gatk_realigner_targets(runner, align_bam,
                                                     ref_file, dbsnp, region,
                                                     out_file, deep_coverage)
        realign_bam = gatk_indel_realignment(runner, align_bam, ref_file,
                                             realign_target_file, region,
                                             out_file, deep_coverage)
        # No longer required in recent GATK (> Feb 2011) -- now done on the fly
        # realign_sort_bam = runner.run_fn("picard_fixmate", realign_bam)
        return realign_bam
    elif out_file:
        shutil.copy(align_bam, out_file)
        return out_file
    else:
        return align_bam

def _has_aligned_reads(align_bam, region):
    """Check if the aligned BAM file has any reads in the region.
    """
    has_items = True
    if region is not None:
        has_items = False
        with closing(pysam.Samfile(align_bam, "rb")) as cur_bam:
            for item in cur_bam.fetch(region):
                has_items = True
                break
    return has_items

# ## High level functionality to run realignments in parallel

def parallel_realign_sample(sample_info, parallel_fn):
    """Realign samples, running in parallel over individual chromosomes.
    """
    to_process = []
    finished = []
    for x in sample_info:
        if x[0]["config"]["algorithm"]["snpcall"]:
            to_process.append(x)
        else:
            finished.append(x)
    if len(to_process) > 0:
        file_key = "work_bam"
        split_fn = split_bam_by_chromosome("-realign.bam", file_key)
        processed = parallel_split_combine(to_process, split_fn, parallel_fn,
                                           "realign_sample", "combine_bam",
                                           file_key, ["config"])
        finished.extend(processed)
    return finished

def realign_sample(data, region=None, out_file=None):
    """Realign sample BAM file at indels.
    """
    log.info("Realigning %s with GATK" % str(data["name"]))
    if data["config"]["algorithm"]["snpcall"]:
        sam_ref = data["sam_ref"]
        config = data["config"]
        data["work_bam"] = gatk_realigner(data["work_bam"], sam_ref, config,
                                          configured_ref_file("dbsnp", config, sam_ref),
                                          region, out_file)
    return [data]
