import sys
import math
import wave
import struct
from scipy.fftpack import fft, fftshift, fftfreq
import array
import bluefile

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

def spectrogram(samprate, data, fftsize=2048, log_scale=True, overlap=0.0, removeDC=False):
    epsilon = 1e-20
    nsamp = len(data)
    first_samp = 0
    last_samp = fftsize
    step = int(fftsize*(1.0-overlap))
    
    samp_period = 1.0 / samprate
    freqs = fftfreq(fftsize, samp_period)
    hdr = bluefile.header(type=2000, format='SF')
    hdr['xstart'] = freqs[0]
    hdr['xdelta'] = freqs[1]-freqs[0]
    hdr['ystart'] = 0.0
    hdr['ydelta'] = step*samp_period
    hdr['subsize'] = fftsize/2
    
    frames = []
    while last_samp < nsamp:
        chunk = data[first_samp:last_samp]
        
        if removeDC:
            dc = float(sum(chunk))/len(chunk)
            chunk = [x - dc for x in chunk]
        
        if len(chunk) < fftsize:
            chunk.extend([0.0]*(fftsize-len(chunk)))
        chunk = fft(chunk)[:fftsize/2]
        
        if log_scale:
            try:
                frame = array.array('f', [10*math.log((x*x.conj()).real+epsilon) for x in chunk])
            except:
                print x, x*x.conj()
                raise
        else:
            frame = array.array('f', [(x*x.conj()).real for x in chunk])
        
        frames.append(frame)

        first_samp += step
        last_samp = first_samp + fftsize
    
    print len(frames), len(frames[0]), len(frames[-1])
    return hdr, frames

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option('-l', '--log', action='store_true', dest='log_scale', default=False, help='Convert to log scale')
    parser.add_option('-f', '--fftsize', action='store', dest='fftsize', default=2048, type='int', help='fft size to use')
    parser.add_option('-o', '--overlap', action='store', dest='overlap', default=0.0, type='float', help='overlap')
    options, args = parser.parse_args()
    if len(args) < 2:
        print 'Must pass input and output filenames on command line'
        sys.exit(-1)
    fname = args[0]
    outname = args[1]
    success, samprate, data = read_simple_wave(fname)
    if not success:
        sys.exit(-1)
    print 'samprate: %s' % samprate
    hdr, frames = spectrogram(samprate, data,
                              log_scale=options.log_scale,
                              fftsize=options.fftsize,
                              overlap=options.overlap)
    bluefile.write(outname, hdr=hdr, data=frames)



