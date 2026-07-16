#!/bin/bash
#SBATCH --job-name=Assembler_thr_HTC040
#SBATCH --partition=HTC040
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --time=1:00:00
#

# !! ADJUST the following line according to your needs: working-subdirectory
cd /nfs/datanumerik/weimann/test/projects/Kaskade7.3/tests/assembler_threading_benchmark
subdir=test10HTC040

mkdir -p $subdir
cd $subdir

# !! ADJUST the following line according to your needs: maximum number of threads
maxthreads=10
export maxthreads
#cmd="srun -B *:*:*  --ntasks=1 ../heat"
cmd="../heat"
rm -f oneThread.time
touch assembsummary.stat
echo "================ order 6 ======================="
for ((i=1;i<=2;i++)); do
  $cmd --init 1 --order 6 --refinements 7 --solver.type iterate --solve false --threads $i 2> /dev/null
done
for ((i=1;i<=$maxthreads;i++)); do
  $cmd --init 0 --order 6 --refinements 7 --solver.type iterate --solve false --threads $i 2> /dev/null
done
mv assembsummary.stat assembsummary_order6.stat

rm -f oneThread.time
touch assembsummary.stat
echo "================ order 5 ======================="
for ((i=1;i<=2;i++)); do
  $cmd --init 1 --order 5 --refinements 7 --solver.type iterate --solve false --threads $i 2> /dev/null
done
for ((i=1;i<=$maxthreads;i++)); do
  $cmd --init 0 --order 5 --refinements 7 --solver.type iterate --solve false --threads $i 2> /dev/null
done
mv assembsummary.stat assembsummary_order5.stat

rm -f oneThread.time
touch assembsummary.stat
echo "================ order 4 ======================="
for ((i=1;i<=2;i++)); do
  $cmd --init 1 --order 4 --refinements 8 --solver.type iterate --solve false --threads $i 2> /dev/null
done
for ((i=1;i<=$maxthreads;i++)); do
  $cmd --init 0 --order 4 --refinements 8 --solver.type iterate --solve false --threads $i 2> /dev/null
done
mv assembsummary.stat assembsummary_order4.stat

rm -f oneThread.time
touch assembsummary.stat
echo "================ order 3 ======================="
for ((i=1;i<=2;i++)); do
  $cmd --init 1 --order 3 --refinements 8 --solver.type iterate --solve false --threads $i 2> /dev/null
done
for ((i=1;i<=$maxthreads;i++)); do
  $cmd --init 0 --order 3 --refinements 8 --solver.type iterate --solve false --threads $i 2> /dev/null
done
mv assembsummary.stat assembsummary_order3.stat

rm -f oneThread.time
touch assembsummary.stat
echo "================ order 2 ======================="
for ((i=1;i<=2;i++)); do
  $cmd --init 1 --order 2 --refinements 9 --solver.type iterate --solve false --threads $i 2> /dev/null
done
for ((i=1;i<=$maxthreads;i++)); do
  $cmd --init 0 --order 2 --refinements 9 --solver.type iterate --solve false --threads $i 2> /dev/null
done
mv assembsummary.stat assembsummary_order2.stat

rm -f oneThread.time
touch assembsummary.stat
echo "================ order 1 ======================="
for ((i=1;i<=2;i++)); do
  $cmd --init 1 --order 1 --refinements 10 --solver.type iterate --solve false --threads $i 2> /dev/null
done
for ((i=1;i<=$maxthreads;i++)); do
  $cmd --init 0 --order 1 --refinements 10 --solver.type iterate --solve false --threads $i 2> /dev/null
done
mv assembsummary.stat assembsummary_order1.stat

/usr/bin/gnuplot ../graph_var_cpus.gnu
/usr/bin/gnuplot ../graph_var_cpus_base2.gnu

