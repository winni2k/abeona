import binascii
import gzip


def is_gz_file(filepath):
    with open(filepath, 'rb') as test_f:
        return binascii.hexlify(test_f.read(2)) == b'1f8b'


def get_maybe_gzipped_file_handle(*args):
    if is_gz_file(args[0]):
        return gzip.open(*args)
    else:
        return open(*args)
