# analyze_dtmf.py

# import all your libraries
import sys
import math
import wave
import struct
from scipy.fftpack import fft, fftshift, fftfreq

# open the audio file
f = wave.open('./phoneCall.wav')
nchannels, sampwidth, fs, nsamp, comptype, compname = f.getparams()
if (nchannels, comptype, compname) != (1, 'NONE', 'not compressed'):
    print 'Failed to read input file.  Params: %s' % f.getparams()
    sys.exit(-1)
if not sampwidth in (1, 2, 4, 8):
    print 'Failed to read input file.  Invalid sampwidth: %s' % sampwidth
    sys.exit(-1)
unpack_type = (None, 'b', 'h', None, 'i', None, None, None, 'q')[sampwidth]
unpack_str = '<%s%s' % (nsamp, unpack_type)

# get the sampling rate (Hz)
fs = f.getframerate()

# compute the sampling period
xdelta = 1.0 / fs

# get number of samples in audio file
nsamp = f.getnframes()

# here we're going to load the data in, one sample at a time. Then we'll compute the mean
# and then we'll subtract the mean off from every point. In this way we'll ensure the 
# signal has zero mean. This isn't strictly necessary - we could just load the data in.

data = struct.unpack(unpack_str, f.readframes(nsamp))
dmean = sum(data)/len(data)     # compute the signal mean
data = [x-dmean for x in data]  # subtract off the mean

# We need to pick an FFT size that will get us an FFT bin size that will allow us
# to distinguish between the two closest tones (697 Hz and 770 Hz), or smaller
# If our FFT bin size (the freq distance between points in the FFT output) is
# larger than about half the minimum distance, it will be difficult to see which
# tone is active on the spectrogram.
max_bin = 25 # Hz
fft_size = fs / max_bin

# FFT math tends to work better when each slice has a length that is a power of 2.
# this chunk of code figures out what the smallest slice size is that matches nslice
# and that is also a power of two. Strictly speaking, this isn't necessary - we could have
# just started by hard-coding a slice size that was a power of 2.
n = int(math.log(fft_size)/math.log(2)+1)
fft_size = int(math.pow(2, n))

# The frame-time or "y-delta" is the time between frames of the spectrogram.  It's best if these
# are short enough so that at least one frame has the tones up for the entire frame (so as to
# avoid the on/off transitions messing with where the energy shows up)
#
ydelta = fft_size / float(fs)
#
# Turns out for the max_bin size we need ydelta is small enough, about 50 ms, so we don't
# have to worry about what to do if this were too large.  The motivated student should
# consider how to deal with a ydelta that is too small.

# this returns the _frequencies_ that go with the FFT. This is one of those pieces of magic
# that's handled for you inside "myFFT.m"
freqs = fftfreq(fft_size, xdelta)

# find the indicies of the "freqs" that correspond to the first and last frequencies of
# the low and high ranges of DTMF freqs. According to wikipedia, the "low" freqs are 697-941Hz
# and the high freqs are 1209-1633Hz. A little extra on either end is always a good idea.

low_start  = None
low_stop   = None
high_start = None
high_stop  = None

for idx, f in enumerate(freqs):
    if low_start is None and f > 600:
        low_start = idx
    elif low_stop is None and f > 1050:
        low_stop = idx
    elif high_start is None and f > 1100:
        high_start = idx
    elif high_stop is None and f > 1800:
        high_stop = idx

def lookup_freqs(f1, f2):
    # this is a subroutine that works out which telephone digit is most likely dialed by freqs
    # f1 and f2
    # note that the default index values of '3' are returned if the frequency doesn't match one of
    # the target values by at least 10000Hz.

    freqs1 = [697, 770, 852, 941]
    freqs2 = [1209, 1336, 1477, 1633]
    dialed = [['1', '2', '3', 'A'], ['4', '5', '6', 'B'], ['7', '8', '9', 'C'], ['*', '0', '#', 'D']]
    
    # find the value of freqs1 that is closest to f1
    min_dist = (10000, 3)
    for i, f in enumerate(freqs1):
        dist = math.fabs(f-f1)
        if dist < min_dist[0]:
            min_dist = (dist, i)
    ind1 = min_dist[1]

    # find the value of freqs2 that is closest to f2    
    min_dist = (10000, 3)
    for i, f in enumerate(freqs2):
        dist = math.fabs(f-f2)
        if dist < min_dist[0]:
            min_dist = (dist, i)
    ind2 = min_dist[1]

    # use ind1 and ind2 to 'lookup' the digit in the 'dialed' array
    return dialed[ind1][ind2]

threshold = 1000000.0

# init first and last samples of the very first slice
first_samp = 0
last_samp = fft_size

# iterate through the slices
while last_samp < nsamp:

    try:
        # grab the data of the current chunk
        chunk = data[first_samp:last_samp]

        # take its fft ; keep only the first half of the samples (eg ignore negative freqs)
        fchunk = fft(chunk)[:fft_size/2]

        # take its magnitude
        mchunk = [math.sqrt((x*x.conj()).real) for x in fchunk]

        csum = 0
        c2sum = 0
        # compute sum and sum-of-squares
        for x in mchunk:
            csum += x
            c2sum += x*x
        # compute the mean magnitude value
        cmean = csum/len(mchunk)
        # compute the magnitude's standard deviation
        cstd = math.sqrt(c2sum/len(mchunk) - cmean*cmean)

        # here we create the "line" of ascii characters to display the FFT
        line = []
        for i in xrange(low_start, high_stop):
            # iterate through the complete range of DTMF freqs
            # if the FFT magnitude is more than 3 standard devs over the
            # mean, plot a *, else plot a -
            if mchunk[i] > cmean+3*cstd:
                line.append('*')
            elif mchunk[i] > cmean+2*cstd:
                line.append('^')
            elif mchunk[i] > cmean+1*cstd:
                line.append('-')
            else:
                line.append(' ')

        # find the maximum FFT magnitude in the "low" DTMF frequency range.
        # store both the magnitude and the corresponding frequency
        low_max = (0, 0)
        for i in xrange(low_start, low_stop):
            if mchunk[i] > low_max[0]:
                low_max = (mchunk[i], freqs[i])

        # find the maximum FFT magnitude in the "high" DTMF frequency range.
        # store both the magnitude and the corresponding frequency
        high_max = (0, 0)
        for i in xrange(high_start, high_stop):
            if mchunk[i] > high_max[0]:
                high_max = (mchunk[i], freqs[i])

        # print the ascii fft along with the sample number of the start of the slice
        # and the low_max and high_max values
        print "%06d" % first_samp, ''.join(line),
        print "low=(%7.2f, %7.2f)" % low_max, "," ,
        print "high=(%7.2f, %7.2f)" % high_max,

        # if both low max and high max magnitudes exceed pre-programmed threshold, 
        # run the digit lookup and print the result
        if low_max[0] > threshold and high_max[0] > threshold:
            print ", digit dialed: ", lookup_freqs(low_max[1], high_max[1])
        else:
            print ''


    finally:
        # iterate the sample indicies so we're all cued up for the next slice
        first_samp += fft_size
        last_samp = first_samp + fft_size


