#!/usr/bin/env python
# Checks the uptime for each node the job is running on
import os
import sys
import subprocess
import re

if len(sys.argv) < 2:
    print >>sys.stderr, "Usage: <job-id>"
    sys.exit(0)

output = subprocess.Popen(["qstat", "-f", sys.argv[1]], stdout=subprocess.PIPE).communicate()[0]
lines = output.split('\n')
hosts = ""
in_hosts = False
for l in lines:
    l = l.strip()
    if l.find("exec_host =") >= 0:
        in_hosts = True
        hosts += l[len("exec_host = "):]
        continue
    if in_hosts and l.find(" = ") >= 0: in_hosts = False
    if in_hosts: hosts += l

hosts = re.sub(r'/[0-9]+', '', hosts)
hosts = sorted(list(set(hosts.split('+'))))
# ssh to all the nodes in parallel
num_kids = 0
for host in hosts:
    pid = os.fork()
    if pid == 0:
        output = subprocess.Popen(["ssh", host, "uptime"], stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0]
        print "%s: %s" % (host, output.strip())
	sys.exit(0)
    else:
        num_kids += 1

while num_kids > 0:
    os.wait()
    num_kids -= 1
