start_time=$(date +%s)

python HB83.py --velocity 4.3 --n_b 100 --sim_per_b 100 --output /home/ybadoux/Documents/MRP/HB83_test2
python HB83.py --velocity 4.98 --n_b 100 --sim_per_b 100 --output /home/ybadoux/Documents/MRP/HB83_test2
python HB83.py --velocity 5.78 --n_b 100 --sim_per_b 100 --output /home/ybadoux/Documents/MRP/HB83_test2
python HB83.py --velocity 6.67 --n_b 100 --sim_per_b 100 --output /home/ybadoux/Documents/MRP/HB83_test2
python HB83.py --velocity 7.71 --n_b 100 --sim_per_b 100 --output /home/ybadoux/Documents/MRP/HB83_test2
python HB83.py --velocity 8.93 --n_b 100 --sim_per_b 100 --output /home/ybadoux/Documents/MRP/HB83_test2
python HB83.py --velocity 10.32 --n_b 100 --sim_per_b 100 --output /home/ybadoux/Documents/MRP/HB83_test2
python HB83.py --velocity 11.95 --n_b 100 --sim_per_b 100 --output /home/ybadoux/Documents/MRP/HB83_test2
python HB83.py --velocity 13.83 --n_b 100 --sim_per_b 100 --output /home/ybadoux/Documents/MRP/HB83_test2
python HB83.py --velocity 16 --n_b 100 --sim_per_b 100 --output /home/ybadoux/Documents/MRP/HB83_test2

end_time=$(date +%s)
runtime=$((end_time - start_time))
echo "Total runtime: $runtime seconds"