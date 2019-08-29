from utils import get_transcript
from read_data import transcripts

from cnv import CNV, CNVRecord

c1 = CNVRecord('11', 2797090, 2869333, 'DEL')
t1 = get_transcript('NM_000218.2', transcripts)
cnv1 = CNV(c1, t1)

c2 = CNVRecord('17', 7572827, 7574133, 'DEL')
t2 = get_transcript('NM_000546.5', transcripts)
cnv2 = CNV(c2, t2)

