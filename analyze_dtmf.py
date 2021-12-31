import math
import wave
import struct

from scipy.fftpack import fft, fftshift, fftfreq

f = wave.open('./phoneCall.wav')
print f.getparams()
xdelta = 1.0 / f.getframerate()
nsamp = f.getnframes()
strsamplen = len(str(nsamp))
strsamp = '%'
strsamp += '%s' % strsamplen
strsamp += 's'
data = []
for i in xrange(nsamp):
    data.append(struct.unpack("<b", f.readframes(1))[0])
dmean = sum(data)/len(data)
data = [x-dmean for x in data]

nslice = 200
slice_size = nsamp/nslice
print slice_size
n = 8
while math.pow(2, n) < slice_size:
    n+=1
slice_size = int(math.pow(2, n))
print slice_size

low_start = None
low_stop = None
high_start = None
high_stop = None

freqs = fftfreq(slice_size, xdelta)
for idx, f in enumerate(freqs):
    if low_start is None and f > 600:
        low_start = idx
    elif low_stop is None and f > 1050:
        low_stop = idx
    elif high_start is None and f > 1100:
        high_start = idx
    elif high_stop is None and f > 1800:
        high_stop = idx

print low_start, freqs[low_start]
print low_stop, freqs[low_stop]
print high_start, freqs[high_start]
print high_stop, freqs[high_stop]

def lookup_freqs(f1, f2):
    freqs1 = [697, 770, 852, 941]
    freqs2 = [1209, 1336, 1477, 1633]
    min_dist = (10000, 3)
    for i, f in enumerate(freqs1):
        dist = math.fabs(f-f1)
        if dist < min_dist[0]:
            min_dist = (dist, i)
    f1 = min_dist[1]
    min_dist = (10000, 3)
    for i, f in enumerate(freqs2):
        dist = math.fabs(f-f2)
        if dist < min_dist[0]:
            min_dist = (dist, i)
    f2 = min_dist[1]
    dialed = [['1', '2', '3', 'A'], ['4', '5', '6', 'B'], ['7', '8', '9', 'C'], ['*', '0', '#', 'D']]
    return dialed[f1][f2]

threshold = 4000.0

first_samp = 0
last_samp = slice_size
while last_samp < nsamp:
    try:
        chunk = data[first_samp:last_samp]
        fchunk = fft(chunk)[:slice_size/2]
        mchunk = [math.sqrt((x*x.conj()).real) for x in fchunk]
        csum = 0
        c2sum = 0
        for x in mchunk:
            csum += x
            c2sum += x*x
        cmean = csum/len(mchunk)
        cstd = math.sqrt(c2sum/len(mchunk) - cmean*cmean)
        line = []
        
        for i in xrange(low_start, high_stop):
            if mchunk[i] > cmean+3*cstd:
                line.append('*')
            elif mchunk[i] > cmean+2*cstd:
                line.append('+')
            elif mchunk[i] > cmean+cstd:
                line.append('-')
            elif mchunk[i] > cmean+0.5*cstd:
                line.append('.')
            else:
                line.append(' ')
        low_max = (0, 0)
        for i in xrange(low_start, low_stop):
            if mchunk[i] > low_max[0]:
                low_max = (mchunk[i], freqs[i])
        high_max = (0, 0)
        for i in xrange(high_start, high_stop):
            if mchunk[i] > high_max[0]:
                high_max = (mchunk[i], freqs[i])
        max_low  = '(%.1f, %.1f) ' % low_max
        max_high = '(%.1f, %.1f)' % high_max
        print strsamp % first_samp, ''.join(line), max_low, max_high,
        if low_max[0] + high_max[0] > 2*threshold:
            print ":", lookup_freqs(low_max[1], high_max[1])
        else:
            print ''
    finally:
        first_samp += slice_size
        last_samp = first_samp + slice_size

# 605-475-6968

