import sys
import math
import wave
import struct
import datetime
from scipy.fftpack import fft, fftshift, fftfreq

def read_simple_wave(fname):
    f = wave.open(fname, 'r')
    nchannels, sampwidth, framerate, nframes, comptype, compname = f.getparams()
    if (nchannels, comptype, compname) != (1, 'NONE', 'not compressed'):
        print 'Failed to read input file.  Params: %s' % f.getparams()
        return False, None, None
    if not sampwidth in (1, 2, 4, 8):
        print 'Failed to read input file. Invalid sampwidth: %s' % sampwidth 
        return False, None, None
    unpack_type = (None, 'b', 'h', None, 'i', None, None, None, 'q')[sampwidth]
    unpack_str = '<%s%s' % (nframes, unpack_type)
    raw_data = f.readframes(nframes)
    data = struct.unpack(unpack_str, raw_data)
    return True, framerate, data

def ascii_spectrogram(samprate, data, fftsize=2048, log_scale=True, cols=128):
    vals = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    nvals = len(vals)
    epsilon = 1e-20

    xdelta = 1.0 / samprate

    nsamp = len(data)
    strsamplen = len(str(nsamp))

    first_samp = 0
    last_samp = fftsize
    bins_per = float(fftsize/2) / cols
    frames = []
    all_vals = []
    while last_samp < nsamp:
        chunk = data[first_samp:last_samp]
        dc = float(sum(chunk))/len(chunk)
        dchunk = [x - dc for x in chunk]
        fchunk = fft(dchunk)[:fftsize/2]
        rchunk = [(x*x.conj()).real for x in fchunk]
        bin1 = 0
        bchunk = []
        for col in xrange(cols):
            bin0 = bin1
            bin1 = int((col+1)*bins_per)
            val = sum(rchunk[bin0:bin1]) / (bin1-bin0)
            bchunk.append(val)
        if log_scale:
            try:
                mchunk = [10*math.log(x+epsilon) for x in bchunk]
            except:
                print x, x*x.conj()
                raise
            frames.append((first_samp, mchunk))
            all_vals.extend(mchunk)
        else:
            frames.append((first_samp, bchunk))
            all_vals.extend(bchunk)
        first_samp += fftsize
        last_samp = first_samp + fftsize
    
    all_vals.sort()
    if log_scale:
        median = all_vals[len(all_vals)/2]
        high_val = all_vals[-1]
        vrange = (high_val-median)/2
        low_val = high_val - vrange
    else:
        median = all_vals[len(all_vals)/2]
        low_val = 100*median
        high_val = all_vals[-1]
        vrange = high_val - low_val
    print median, low_val, high_val, vrange
    
    for first_samp, chunk in frames:
        line = []
        for val in chunk:
            if val > low_val:
                letter = min(35, max(0, int((val-low_val)/vrange*nvals)))
                line.append(vals[letter])
            else:
                line.append(' ')
        tstamp = '%-11s' % str(datetime.timedelta(seconds=(first_samp*xdelta)))[:11]
        print tstamp, '|', ''.join(line), '|'
                
        

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option('-l', '--log', action='store_true', dest='log_scale', default=False, help='Convert to log scale')
    parser.add_option('-f', '--fftsize', action='store', dest='fftsize', default=2048, type='int', help='fft size to use')
    parser.add_option('-c', '--cols', action='store', dest='cols', default=128, type='int', help='display columns')
    options, args = parser.parse_args()
    if len(args) < 1:
        print 'Must pass input filename on command line'
        sys.exit(-1)
    fname = args[0]
    success, framerate, data = read_simple_wave(fname)
    if not success:
        sys.exit(-1)
    print 'framerate: %s' % framerate
    ascii_spectrogram(framerate, data, log_scale=options.log_scale, fftsize=options.fftsize, cols=options.cols)


