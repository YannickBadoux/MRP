start_time=$(date +%s)

python acs_binary.py --velocity 1.02 --n_sim 1000 --output /home/ybadoux/Documents/MRP/acs_1000
python acs_binary.py --velocity 4.3 --n_sim 1000 --output /home/ybadoux/Documents/MRP/acs_1000
python acs_binary.py --velocity 2.07 --n_sim 1000 --output /home/ybadoux/Documents/MRP/acs_1000
python acs_binary.py --velocity 5.78 --n_sim 1000 --output /home/ybadoux/Documents/MRP/acs_1000
python acs_binary.py --velocity 1.79 --n_sim 1000 --output /home/ybadoux/Documents/MRP/acs_1000
python acs_binary.py --velocity 6.67 --n_sim 1000 --output /home/ybadoux/Documents/MRP/acs_1000
python acs_binary.py --velocity 2.40 --n_sim 1000 --output /home/ybadoux/Documents/MRP/acs_1000
python acs_binary.py --velocity 10.32 --n_sim 1000 --output /home/ybadoux/Documents/MRP/acs_1000
python acs_binary.py --velocity 3.21 --n_sim 1000 --output /home/ybadoux/Documents/MRP/acs_1000
python acs_binary.py --velocity 13.83 --n_sim 1000 --output /home/ybadoux/Documents/MRP/acs_1000

end_time=$(date +%s)
runtime=$((end_time - start_time))
echo "Total runtime: $runtime seconds"