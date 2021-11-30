MATSIZE=$1

for KVAL in {1..10}
do
echo ./build/krylov.exe -t import -n $MATSIZE -i $MATSIZE -r $MATSIZE -k $KVAL -o 1
./build/krylov.exe -t import -n $MATSIZE -i $MATSIZE -r $MATSIZE -k $KVAL -o 1
gprof ./build/krylov.exe gmon.out > ./results/n${MATSIZE}/output_nosimd${KVAL}.txt
echo ""
done