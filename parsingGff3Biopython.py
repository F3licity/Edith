import pprint
from BCBio.GFF import GFFExaminer
 
in_file = "Nagalakshmi_2008_UTRs.gff3"
examiner = GFFExaminer()
in_handle = open(in_file)
pprint.pprint(examiner.parent_child_map(in_handle))
in_handle.close()

from BCBio import GFF
 
in_file = "Nagalakshmi_2008_UTRs.gff3"
 
in_handle = open(in_file)
for rec in GFF.parse(in_handle):
    print rec
in_handle.close()
