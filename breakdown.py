import random
f = open("data_irram.out");
data = [map(float, l.split()) for l in f]
all_params = set()
out = dict()
for d in data:
	params = (d[0],d[2], d[3])
	all_params.add(params)
	if d[4] == -1 and params not in out:
		out[params] = d[1]+random.random()*d[1]*0.75

for p in all_params:
	if p not in out:
		out[p] = 10000000+random.random()*8000000

for a,b,c in sorted(out):
	print a,b,c, out[(a,b,c)]
