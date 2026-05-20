#! /usr/bin/env bash

docker run --mount type=bind,src=./output_data,dst=/app/output_data nbarron00/nn-grid-search-general ./training_data/data_bias_e5_train.csv 'FValue' ./output_data/dummy_output 42 20 8 100000 200 ./training_data/data_bias_e5_val.csv ./output_data/test_model 1 0.001