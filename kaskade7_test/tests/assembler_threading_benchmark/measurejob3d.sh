#!/bin/bash
#SBATCH --job-name=Assembler_thr_HTC050_3d
#SBATCH --partition=HTC050
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --time=72:00:00
#

# !! ADJUST the following line according to your needs: working-subdirectory
cd /nfs/datanumerik/weimann/test/projects/Kaskade7.3/tests/assembler_threading_benchmark
subdir=test3d48HTC050

mkdir -p $subdir
cd $subdir

# !! ADJUST the following line according to your needs: maximum number of threads
maxthreads=48
export maxthreads
#cmd="srun -B *:*:*  --ntasks=1 ../heat3d"
cmd="../heat3d"

rm -f oneThread3d.time
touch assembsummary3d.stat
echo "================ order 6 ======================="
for ((i=1;i<=2;i++)); do
  $cmd --init 1 --order 6 --refinements 3 --heapsize 8192 --solver.type iterate --solve false --threads $i 2> /dev/null
done
for ((i=1;i<=$maxthreads;i++)); do
  $cmd --init 0 --order 6 --refinements 3 --heapsize 8192 --solver.type iterate --solve false --threads $i 2> /dev/null
done
mv assembsummary3d.stat assembsummary3d_order6.stat

rm -f oneThread3d.time
touch assembsummary3d.stat
echo "================ order 5 ======================="
for ((i=1;i<=2;i++)); do
  $cmd --init 1 --order 5 --refinements 4 --heapsize 8192 --solver.type iterate --solve false --threads $i 2> /dev/null
done
for ((i=1;i<=$maxthreads;i++)); do
  $cmd --init 0 --order 5 --refinements 4 --heapsize 8192 --solver.type iterate --solve false --threads $i 2> /dev/null
done
mv assembsummary3d.stat assembsummary3d_order5.stat

rm -f oneThread3d.time
touch assembsummary3d.stat
echo "================ order 4 ======================="
for ((i=1;i<=2;i++)); do
  $cmd --init 1 --order 4 --refinements 4 --heapsize 8192 --solver.type iterate --solve false --threads $i 2> /dev/null
done
for ((i=1;i<=$maxthreads;i++)); do
  $cmd --init 0 --order 4 --refinements 4 --heapsize 8192 --solver.type iterate --solve false --threads $i 2> /dev/null
done
mv assembsummary3d.stat assembsummary3d_order4.stat

rm -f oneThread3d.time
touch assembsummary3d.stat
echo "================ order 3 ======================="
for ((i=1;i<=2;i++)); do
  $cmd --init 1 --order 3 --refinements 4 --heapsize 8192 --solver.type iterate --solve false --threads $i 2> /dev/null
done
for ((i=1;i<=$maxthreads;i++)); do
  $cmd --init 0 --order 3 --refinements 4 --heapsize 8192 --solver.type iterate --solve false --threads $i 2> /dev/null
done
mv assembsummary3d.stat assembsummary3d_order3.stat

rm -f oneThread3d.time
touch assembsummary3d.stat
echo "================ order 2 ======================="
for ((i=1;i<=2;i++)); do
  $cmd --init 1 --order 2 --refinements 5 --heapsize 8192 --solver.type iterate --solve false --threads $i 2> /dev/null
done
for ((i=1;i<=$maxthreads;i++)); do
  $cmd --init 0 --order 2 --refinements 5 --heapsize 8192 --solver.type iterate --solve false --threads $i 2> /dev/null
done
mv assembsummary3d.stat assembsummary3d_order2.stat

rm -f oneThread3d.time
touch assembsummary3d.stat
echo "================ order 1 ======================="
for ((i=1;i<=2;i++)); do
  $cmd --init 1 --order 1 --refinements 6 --heapsize 8192 --solver.type iterate --solve false --threads $i 2> /dev/null
done
for ((i=1;i<=$maxthreads;i++)); do
  $cmd --init 0 --order 1 --refinements 6 --heapsize 8192 --solver.type iterate --solve false --threads $i 2> /dev/null
done
mv assembsummary3d.stat assembsummary3d_order1.stat

/usr/bin/gnuplot ../graph_3d_var_cpus.gnu
/usr/bin/gnuplot ../graph_3d_var_cpus_base2.gnu

