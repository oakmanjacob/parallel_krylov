MATSIZE=$1

for KVAL in {1..10}
do
echo ./build/krylov.exe -m import -n $MATSIZE -i $MATSIZE -r 24 -k $KVAL -o 2
./build/krylov.exe -m import -n $MATSIZE -i $MATSIZE -r 24 -k $KVAL -o 2
gprof ./build/krylov.exe gmon.out > ./results/n${MATSIZE}/output_nosimd${KVAL}.txt
echo ""

echo ./build/krylov.exe -m import -n $MATSIZE -i $MATSIZE -r 24 -k $KVAL -o 4
./build/krylov.exe -m import -n $MATSIZE -i $MATSIZE -r 24 -k $KVAL -o 4
gprof ./build/krylov.exe gmon.out > ./results/n${MATSIZE}/output_nosimd${KVAL}.txt
echo ""
done