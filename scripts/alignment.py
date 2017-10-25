class Alignment:
    def __init__(self, qry, ref, qrystartpos, qryendpos, refstartpos, refendpos, orientation, confidence, hitenum, qrylen, reflen, channel, alignmentstring): 
        self.qry = qry
        self.ref = ref
        self.qrystartpos = qrystartpos  
        self.qryendpos = qryendpos
        self.refstartpos = refstartpos
        self.refendpos = refendpos
        self.orientation = orientation
        self.confidence = confidence
        self.qrylen = qrylen
        self.reflen = reflen
        self.hitenum = hitenum
        self.channel = channel
        self.alignmentstring = alignmentstring
        self.start = 0 # will be filled later
        self.end = 0   # will be filled later
        self.qry_left_overlen = 0 # will be filled later
        self.qrt_right_overlen = 0 # will be filled later
        self.chimeric = False # will be change later
        self.contained = False # will be change later

    def __str__(self):
        return 'ref %d refstart %.1f refend %.1f qry %d qrylen %.1f qrystart %.1f qryend %.1f %s conf %.2f start %.2f end %2.f qry_left_overlen %2.f qrt_right_overlen %2.f chimeric %s contained %s' % (self.ref, self.refstartpos, self.refendpos, self.qry, self.qrylen, self.qrystartpos, self.qryendpos, self.orientation, self.confidence, self.start, self.end, self.qry_left_overlen, self.qrt_right_overlen, self.chimeric, self.contained)
        #return 'ref %d, reflen %.1f, refstartpos %.1f, refendpos %.1f, qry %d, qrylen %.1f, qrystartpos %.1f, qryendpos %.1f, orientation %s, confidence %.2f' % (self.ref, self.reflen, self.refstartpos, self.refendpos, self.qry, self.qrylen, self.qrystartpos, self.qryendpos, self.orientation, self.confidence)

    def unpack(self):
        return [self.qry, self.ref, self.qrystartpos, self.qryendpos, self.refstartpos, self.refendpos, self.orientation, self.confidence, self.hitenum, self.qrylen, self.reflen, self.channel, self.alignmentstring]



