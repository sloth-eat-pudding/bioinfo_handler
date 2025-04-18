from handler.utils import *
from handler.config_parser import get_config_column

class CommandRunner:
    def __init__(self, ):
        self.commands = []
        server_name = socket.gethostname()
        longphase = None
        ref_path = None
        self.threads = 24
        data_root_path = get_config_column("data_root_path")
        run_root_path = get_config_column("run_root_path")

        if "gpu" in server_name:
            ref_path = f"{data_root_path}/GRCh38_no_alt_analysis_set.fasta"
            self.threads = 24
        elif "bip8" in server_name:
            ref_path = "/mnt/ramdisk/GRCh38_no_alt_analysis_set.fasta"
            # ref_path = "/big8_disk/giab_lsk114_2022.12/alignment/sup/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"
            self.threads = 50
        longphase = f"{run_root_path}/longphase-to"
        pon_path = f"{data_root_path}/PON"

        self.longphase = longphase
        self.ref_path = ref_path
        self.pon_path = pon_path

    def mkdir_folder(self, output_path):
        output_dir = os.path.dirname(output_path)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

    def add_make_command(self, path:str = None):
        if path is None:
            path = self.longphase
        self.commands.insert(0, ["make clean -C " + path])
        self.commands.insert(1, ["make -j 20 -C " + path])
        return self
    
    def add_phase_command(self, bam_path, vcf_path, output_path, indel=False, background=False, test=False, tmp_num=""):
        self.mkdir_folder(output_path)
        command = [
            # "valgrind --leak-check=full --show-leak-kinds=all -v -s ",
            f"/bip8_disk/zhenyu112/longphase{tmp_num}",
            # f"{self.longphase}/longphase",
            "phase",
            "-b", bam_path,
            "-r", self.ref_path,
            "-s", vcf_path,
            "-t", str(self.threads),
            "--ont",
            "-o", output_path,
            "--ge",
            "--lge",
            "--sge",
            "--loh",
            "--cnv",
            # "--somaticConnectAdjacent", "6",
            "--somaticConnectAdjacent", "8",
            # "--pon-file", f"{self.pon_path}/clairs-to_databases/1000g-pon.sites.vcf,{self.pon_path}/clairs-to_databases/CoLoRSdb.GRCh38.v1.1.0.deepvariant.glnexus.af-ge-0.001.vcf,{self.pon_path}/clairs-to_databases/dbsnp.b138.non-somatic.sites.vcf,{self.pon_path}/clairs-to_databases/gnomad.r2.1.af-ge-0.001.sites.vcf",
            "--pon-file", f"{self.pon_path}/clairs-to_databases/1000g-pon.sites.vcf,{self.pon_path}/clairs-to_databases/CoLoRSdb.GRCh38.v1.1.0.deepvariant.glnexus.af-ge-0.001.vcf",
            "--strict-pon-file", f"{self.pon_path}/clairs-to_databases/dbsnp.b138.non-somatic.sites.vcf,{self.pon_path}/clairs-to_databases/gnomad.r2.1.af-ge-0.001.sites.vcf",
            # "--pon-file", f"{self.pon_path}/deepsomatic-to_databases/PON_dbsnp138_gnomad_PB1000g_pon.vcf,{self.pon_path}/deepsomatic-to_databases/AF_pacbio_PON_CoLoRSdb.GRCh38.AF0.05.vcf",
            # "-a", "35"
        ]
        if indel:
            command.append("--indels")
        if background:
            command.append(" &")
        if test:
            command.insert(0, "/usr/bin/time -v ")
        self.commands.append(command)
        return self
    
    def add_tag_command(self, bam_path, vcf_path, output_path, background=False, test=False):
        self.mkdir_folder(output_path)
        command = [
            f"{self.longphase}/longphase",
            "haplotag",
            "-b", bam_path,
            "-r", self.ref_path,
            "-s", vcf_path+".vcf",
            "-t", str(self.threads),
            "-q", "0",
            "-o", output_path,
            # "--region", "chr2"
        ]
        if background:
            command.append(" &")
        if test:
            command.insert(0, "/usr/bin/time -v ")
        self.commands.append(command)
        return self
    
    def add_index_command(self, bam_path, background=False, test=False):
        command = [
            "samtools index -@ 50 ",
            bam_path,
        ]
        if background:
            command.append(" &")
        if test:
            command.insert(0, "/usr/bin/time -v ")
        self.commands.append(command)
        return self
    
    def add_vcf_index_command(self, vcf_path, background=False, test=False):
        command = [
            # f"bgzip -k {vcf_path}.vcf && tabix -p vcf {vcf_path}.vcf.gz",
            f"rm -f {vcf_path}.vcf.gz {vcf_path}.vcf.gz.tbi && bgzip -k {vcf_path}.vcf && tabix -p vcf {vcf_path}.vcf.gz"
        ]
        if background:
            command.append(" &")
        if test:
            command.insert(0, "/usr/bin/time -v ")
        self.commands.append(command)
        return self
    
    def add_merge_command(self, input_chr_file, parser_file, run_mode = "", background=False):
        command = []
        if run_mode == "cat":
            command = [
                "cat", f"{input_chr_file}chr*.{parser_file}", ">", f"{input_chr_file}_{parser_file}"
            ]
        elif run_mode == "awk":
            command = [
                # "bash", "-c",
                "awk",
                # "'$7 != \"0\" || $9 != \"0\" || $11 != \"0\" || $13 != \"0\" || $15 != \"0\" || $17 != \"0\" || $19 != \"0\" || $21 != \"0\" {print $1, $3, $5}'",
                "'$5 != \"-1\" && ($7 != \"0\" || $9 != \"0\" || $11 != \"0\" || $13 != \"0\" || $15 != \"0\" || $17 != \"0\" || $19 != \"0\" || $21 != \"0\") {print $1, $3, $5}'",
                # "'$5 != \"-1\" || $7 != \"0\" || $9 != \"0\" || $11 != \"0\" || $13 != \"0\" || $15 != \"0\" || $17 != \"0\" || $19 != \"0\" || $21 != \"0\" {print $1, $2, $3, $4, $5, $7, $9, $11, $13, $15, $17, $19, $21}'",
                f"{input_chr_file}chr*.{parser_file}", ">", f"{input_chr_file}_{parser_file}"
            ]
        if background:
            command.append(" &")
        self.commands.append(command)
        return self

    def add_command(self, command:list, output_path=None, background=False, test=False):
        if output_path is not None:
            self.mkdir_folder(output_path)
        if background:
            if isinstance(command, list):
                # self.commands.append(command)
                command.append(" &")
            elif isinstance(command, str):
                command += " &"
        if test:
            command.insert(0, "/usr/bin/time -v")
        self.commands.append([command])
        return self

    def run(self, log_file_path:str = None, hide_log=False, check_return=True):
        """
        run and write log
        """
        if log_file_path is None:
            log_file_path = "execution.log"
        if not os.path.exists(log_file_path):
            os.makedirs(os.path.dirname(log_file_path), exist_ok=True)
        # process log file name
        base_name, ext = os.path.splitext(log_file_path)
        counter = 1
        base_name, counter = re.subn(r'(_\d+)?$', '', base_name)
        while os.path.exists(log_file_path):
            log_file_path = f"{base_name}_{counter}{ext}"
            counter += 1

        with open(log_file_path, "w") as log_file:
            for command in self.commands:
                command_str = ' '.join(command)
                logging.info(command_str)
                log_file.write(f"{command_str}\n")
                if not hide_log:
                    process = subprocess.run(command_str, shell=True, check=check_return)
                else:
                    process = subprocess.run(command_str, shell=True, check=check_return, stdout=log_file, stderr=log_file)
                if process.returncode != 0:
                    log_file.write(f"Command '{command}' failed with return code {process.returncode}\n")
        logging.info("run command end")
        self.commands = []
        return self

