"""Evaluate saved mapping CSVs through BReMAP, without research-script imports."""
import argparse
import json
import bremap as br

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--teacher',required=True,help='Teacher cells.csv')
p.add_argument('--student',required=True,help='Student cells.csv')
p.add_argument('--truth',required=True,help='CSV: cell_id, donor_id, label')
p.add_argument('--states',nargs='+',required=True)
p.add_argument('--agreement',type=float,default=.85)
p.add_argument('--minimum-f1',type=float,required=True)
p.add_argument('--minimum-coverage',type=float,required=True)
args=p.parse_args()
report=br.evaluate_student(args.teacher,args.student,truth=args.truth,states=args.states,
    agreement_threshold=args.agreement,maximum_coverage_change=.05,maximum_f1_loss=.05,
    minimum_task_f1=args.minimum_f1,minimum_coverage=args.minimum_coverage)
print(json.dumps(report,indent=2))
