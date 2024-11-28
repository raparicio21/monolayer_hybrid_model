cp inputs/* .  
echo "Random info I want to print"    

mkdir results
/usr/local/matlab/R2021a/bin/matlab -nosplash -nodesktop -r "run('topology.m'); exit;"
tar -czf results.tar.gz results