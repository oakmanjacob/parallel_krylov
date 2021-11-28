MATSIZE=$1

for KVAL in {1..10}
do
echo ./build/krylov.exe -t import -n $MATSIZE -i $MATSIZE -r $MATSIZE -k $KVAL -o
./build/krylov.exe -t import -n $MATSIZE -i $MATSIZE -r $MATSIZE -k $KVAL -o
gprof ./build/krylov.exe gmon.out > ./data/n${MATSIZE}/output${KVAL}.txt
echo ""
done