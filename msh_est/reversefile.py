"""
    reverse the lines of a text file

    should work regardless of file size and ram

    JHey 6/5/2018
"""

import sys
import os


def reverse_readline(filename, buf_size=8192):
    """a generator that returns the lines of a file in reverse order
        found at  https://stackoverflow.com/questions/2301789/read-a-file-in-reverse-order-using-python/26747854
    """
    with open(filename,'rb') as fh:      ## use binary mode,  apparently more accurate when using tell() seek()
        segment = None
        offset = 0
        fh.seek(0, os.SEEK_END)
        file_size = remaining_size = fh.tell()
        while remaining_size > 0:
            offset = min(file_size, offset + buf_size)
            fh.seek(file_size - offset)
            buffer = fh.read(min(remaining_size, buf_size))
##            remaining_size -= buf_size  replace with len(buffer)
            remaining_size -= len(buffer)
            lines = buffer.decode('ascii').split('\n')    ## decode the binary back to ascii
##            lines = buffer.split('\n')
            # the first line of the buffer is probably not a complete line so
            # we'll save it and append it to the last line of the next buffer
            # we read
            if segment is not None:
                # if the previous chunk starts right from the beginning of line
                # do not concact the segment to the last line of new chunk
                # instead, yield the segment first
                if buffer[-1] is not '\n':
                    lines[-1] += segment
                else:
                    yield segment
            segment = lines[0]
            for index in range(len(lines) - 1, 0, -1):
                if len(lines[index]):
                    yield lines[index]
        # Don't yield None if the file was empty
        if segment is not None:
            yield segment


def linecount(filename):
    """
        https://stackoverflow.com/questions/845058/how-to-get-line-count-cheaply-in-python
    """
    f = open(filename, 'rb')
    lines = 0
    buf_size = 1024 * 1024
    read_f = f.raw.read
    buf = read_f(buf_size)
    while buf:
        lines += buf.count(b'\n')
        buf = read_f(buf_size)

    return lines


def reverse_file(fname,frevname,doprint = None):
    if doprint:
        print(fname,"has",linecount(fname),"lines")
    frout = open(frevname,"w")
    i = 0
    lineintervaloutput = 1000000
    for rline in reverse_readline(fname):
        frout.write(rline.strip()+"\n")   # strip in case leftover \n or \r
        i += 1
        if doprint:
            if i == lineintervaloutput * (i//lineintervaloutput):
                print("i","lines written")
    frout.close()

if __name__ == "__main__":
    if len(sys.argv) == 3:
        fname = sys.argv[1]
        frevname = sys.argv[2]
        reverse(fname,frevname,doprint = True)
    else:
        print("usage: python reversefile.py filename1  filname2\n  filename1 is the name of the file to be reversed\n  filename2 is the name of the reversed file")


