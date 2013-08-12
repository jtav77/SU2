# Read lines from file
f = open('cmesh.su2','r')
lines = f.readlines()
f.close()

# Write lines back to file, stripping away leading spaces
f = open('cmesh.su2','w')
for line in lines:
  line2 = line.strip()
  f.write('%s\n'%line2)
f.close()
