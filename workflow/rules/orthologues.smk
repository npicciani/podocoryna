rule orthofinder:
    input:
        directory("resources/orthofinder")
    output:
        directory("results/orthofinder")
    conda:
        "../envs/orthofinder.yaml"
    log:
        "logs/orthofinder/orthofinder.log"
    threads: 20
    shell:
        "orthofinder -f {input} -t {threads} -o {output}"
