MATSIZE=$1

for KVAL in {1..10}
do
echo ./build/krylov.exe -m import -n $MATSIZE -i $MATSIZE -r 100 -k $KVAL -o $2 -t 8
./build/krylov.exe -m import -n $MATSIZE -i $MATSIZE -r 100 -k $KVAL -o $2 -t 8
gprof ./build/krylov.exe gmon.out > ./results/n${MATSIZE}/output_nosimd${KVAL}.txt
echo ""
done