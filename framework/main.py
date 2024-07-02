#!/bin/python3

# main.py
import importlib
import subprocess

def run_simulation(config):
    subprocess.run(["bash", "scripts/run_simulation.sh", config["input_file"], config["output_directory"]])

def main():
    configs = ["configs.config1"]

    for config_name in configs:
        config = importlib.import_module(config_name).config

        print(config)
        
        # Run the simulation
        run_simulation(config)

if __name__ == "__main__":
    main()
