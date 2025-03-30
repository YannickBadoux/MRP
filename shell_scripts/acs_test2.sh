start_time=$(date +%s)

python acs_binary.py --velocity 1.16 --n_sim 100 --output /home/ybadoux/Documents/MRP/acs_100
python acs_binary.py --velocity 4.98 --n_sim 100 --output /home/ybadoux/Documents/MRP/acs_100
python acs_binary.py --velocity 1.55 --n_sim 100 --output /home/ybadoux/Documents/MRP/acs_100
python acs_binary.py --velocity 1.34 --n_sim 100 --output /home/ybadoux/Documents/MRP/acs_100
python acs_binary.py --velocity 7.71 --n_sim 100 --output /home/ybadoux/Documents/MRP/acs_100
python acs_binary.py --velocity 8.93 --n_sim 100 --output /home/ybadoux/Documents/MRP/acs_100
python acs_binary.py --velocity 2.78 --n_sim 100 --output /home/ybadoux/Documents/MRP/acs_100
python acs_binary.py --velocity 11.95 --n_sim 100 --output /home/ybadoux/Documents/MRP/acs_100
python acs_binary.py --velocity 3.72 --n_sim 100 --output /home/ybadoux/Documents/MRP/acs_100
python acs_binary.py --velocity 16.0 --n_sim 100 --output /home/ybadoux/Documents/MRP/acs_100

end_time=$(date +%s)
runtime=$((end_time - start_time))
echo "Total runtime: $runtime seconds"