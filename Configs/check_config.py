import yaml
yaml_file = "/public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err/Pipe/Configs/Config.yaml"
with open(yaml_file, "r") as f:
    config = yaml.safe_load(f.read())
    # print(config)
    ls = (config["denovo_asm"].get("flye", []))
    print(ls)
    print(" ".join(ls))
    pass