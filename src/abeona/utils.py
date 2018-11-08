def is_gz_file(filepath):
    import binascii
    with open(filepath, 'rb') as test_f:
        return binascii.hexlify(test_f.read(2)) == b'1f8b'


def get_maybe_gzipped_file_handle(*args):
    import gzip
    if is_gz_file(args[0]):
        return gzip.open(*args)
    else:
        return open(*args)


def get_fastx_title_seq_generator(fh):
    """Requires a seekable file handle. Needs file to start at beginning of file handle."""
    from Bio.SeqIO.QualityIO import FastqGeneralIterator
    from Bio.SeqIO.FastaIO import SimpleFastaParser

    # determine if the file looks like a fastq
    fh.seek(0)
    is_fastq = True
    try:
        next(FastqGeneralIterator(fh), None)
    except ValueError as e:
        if e.args[0] == "Records in Fastq files should start with '@' character":
            is_fastq = False
        else:
            raise
    fh.seek(0)
    if is_fastq:
        return ((title, seq) for title, seq, qual in FastqGeneralIterator(fh))
    else:
        return SimpleFastaParser(fh)
