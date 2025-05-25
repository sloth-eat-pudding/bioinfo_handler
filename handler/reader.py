from handler.utils import *
from handler.config_parser import get_config_column


@pd.api.extensions.register_dataframe_accessor("custom")
class CustomAccessor:
    def __init__(self, pandas_obj):
        self._obj = pandas_obj

    def sort_chr_pos(self, dictionary: bool = False, inplace: bool = False) -> pd.DataFrame:
        """
        按照染色体和位置排序 DataFrame。

        参数:
        - dictionary (bool): 决定按哪个染色体列排序。
                              如果为 True，按 'CHROM' 排序；否则，按 'CHROM_num' 排序。
        - inplace (bool): 是否在原 DataFrame 上进行排序。默认为 False，返回新的 DataFrame。

        返回:
        - pd.DataFrame: 排序后的 DataFrame，如果 inplace=True，则返回 None。
        """
        # 确保 DataFrame 拥有必要的列
        if 'CHROM' not in self._obj.columns:
            raise KeyError("DataFrame 必须包含 'CHROM' 列。")

        sort_pos_col = "POS1" if "POS1" in self._obj.columns else "POS"

        if sort_pos_col not in self._obj.columns:
            raise KeyError(f"DataFrame 必须包含 '{sort_pos_col}' 列。")

        # 创建排序后的 DataFrame 副本
        df_sorted = self._obj.copy()

        if not dictionary:
            # 使用 .loc 进行赋值，避免 SettingWithCopyWarning
            df_sorted.loc[:, "CHROM_num"] = df_sorted["CHROM"].map(chr_to_num)

        # 执行排序
        sort_chr_col = "CHROM" if dictionary else "CHROM_num"
        df_sorted = df_sorted.sort_values(
            by=[sort_chr_col, sort_pos_col],
            ascending=[True, True]
        )

        if not dictionary:
            # 删除临时列
            df_sorted = df_sorted.drop(columns=["CHROM_num"])

        # 重置索引
        df_sorted = df_sorted.reset_index(drop=True)

        if inplace:
            # 如果是原地排序，修改原 DataFrame
            # self._obj.drop(self._obj.index, inplace=True)
            # self._obj[:] = df_sorted.values
            # return None
            # for column in df_sorted.columns:
            #     self._obj[column] = df_sorted[column].values
            # return None
            for column in df_sorted.columns:
                self._obj.loc[:, column] = df_sorted[column].values
            # 如果 DataFrame 的列有变化（例如新增或删除列），需要处理
            # 删除不在 df_sorted 中的列
            columns_to_drop = set(self._obj.columns) - set(df_sorted.columns)
            if columns_to_drop:
                self._obj.drop(columns=list(columns_to_drop), inplace=True)
            return None
        else:
            return df_sorted

    def rename_merge_and_map(self, mapping, column_name, inplace: bool = True) -> pd.DataFrame:
        """
        将 '_merge' 列的值映射到新的列名，并重命名 '_merge' 列。

        参数:
        - mapping (dict): 映射字典，用于将 '_merge' 列的值替换为其他值。
        - column_name (str): 新列的名称。
        - inplace (bool): 是否原地修改 DataFrame。默认为 True。

        返回:
        - pd.DataFrame: 修改后的 DataFrame。如果 inplace=True，则返回 None。
        """
        if '_merge' not in self._obj.columns:
            raise KeyError("DataFrame 必须包含 '_merge' 列。")

        df_modified = self._obj.copy()
        df_modified['_merge'] = df_modified['_merge'].map(mapping)
        df_modified.rename(columns={'_merge': column_name}, inplace=True)

        if inplace:
            for column in df_modified.columns:
                self._obj.loc[:, column] = df_modified[column].values

            columns_to_drop = set(self._obj.columns) - set(df_modified.columns)
            if columns_to_drop:
                self._obj.drop(columns=list(columns_to_drop), inplace=True)
            return None
        else:
            return df_modified


class BioinfoFileReader:
    """
    - reader(replace_sample = 'COLO-829BL-NovaSeq')
    - filter_by_name(name = 'NonSomatic')
    - select_by_name(name = 'PASS')
    - cut_filter
    - cut_info
    - cut_format_sample(only_save_columns = ['GT'])
    - check_variant_type
    - sort_chr_pos
    """

    def __init__(self, file_path: str = None, columns: list = None, second_column: int = 0):
        # logging.info("BioinfoFileReader---------------------------------------------------\n"
        #              f"\t\t\t       {file_path}")
        logging.info(f"BioinfoFileReader:{file_path}")
        self.file_pattern = file_path

        self.sep = None
        self.columns = columns
        self.dataformat = None
        self.df = None
        self.alt_column_names = None
        self.alt_len_column_names = None
        self.multiple_files = False
        self.file_format_set(second_column)

    def file_format_set(self, second_column: int = 0):
        _, ext = os.path.splitext(self.file_pattern)
        ext = ext.lower().strip(".")
        columns = None
        if ext == "vcf" or ".vcf.gz" in self.file_pattern:
            sep = r"\s+"
            self.dataformat = "vcf"
            # columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT',
            #            'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE']
        elif ext == "bed":
            sep = r"\s+"
            if second_column == 0:
                columns = ["CHROM", "POS1", "POS2", "TAG"]
            elif second_column == 1:
                columns = ["CHROM", "POS1", "POS2", "TAG", "name1", "name2"]
            elif second_column == 2:
                columns = ["CHROM", "POS1", "POS2", "NAME", "SCORE", "STRAND", "THICKSTART", "THICKEND", "ITEM_RGB"]
        elif ext == "block":
            sep = r"\s+"
            columns = ["CHROM", "POS1", "LEN", "POS2"]
        elif ext == "mpileup" or ext == "pileup":
            ext = "pileup"
            sep = r"\s+"
            columns = ["CHROM", "POS", "REF",
                       "COVERAGE", "READ_BASE", "BASE_QUAL"]
        elif ext == "sam":
            sep = r"\t"
        else:
            # sep = r'\t'
            sep = r"\s+"
            if self.columns is None:
                logging.warning(f"格式沒有通用支援")
        if "*" in self.file_pattern:
            self.multiple_files = True
        self.ext = ext
        self.sep = sep
        if self.columns is None:
            self.columns = columns

    def reader(self, replace_sample: str = None):
        logging.debug(f"reader start")
        if not self.multiple_files:
            if self.dataformat == "vcf":
                last_header_line = None
                with open(self.file_pattern, "r") as f:
                    for line in f:
                        # 交界處會是標頭
                        if line.startswith("#"):
                            last_header_line = line
                        else:
                            break
                if last_header_line:
                    self.columns = last_header_line.strip().lstrip("#").split("\t")
                else:
                    logging.warning("last_header_line 為空，給定預設10欄位(須根據實際情況更改)")
                    self.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE']
                    # self.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
                if len(self.columns) != 10:
                    logging.warning(f"vcf columns != 10")
                    logging.debug(f"columns: {self.columns}")
                # 如果沒有SAMPLE，則會用replace_sample或Sample取代，因為有些檔案會是小寫
                if 'SAMPLE' not in self.columns and (replace_sample or 'Sample' in self.columns):
                    self.columns = ["SAMPLE" if col in {replace_sample, 'Sample'} else col for col in self.columns]
                logging.debug(f"columns: {self.columns}")
            self.df = pd.read_csv(self.file_pattern, sep=self.sep, comment="#", names=self.columns)
        else:
            file_paths = glob.glob(self.file_pattern)
            for file_path in file_paths:
                df = pd.read_csv(file_path, sep=self.sep, comment="#", names=self.columns)
                self.df = pd.concat([self.df, df], ignore_index=True)
        logging.debug(f"reader end")
        return self

    def filter_by_name(self, name: str):
        logging.debug(f"filter_by_columns start {name}")
        original_count = len(self.df)
        self.df = self.df[~self.df["FILTER"].str.contains(name, na=False)]
        filtered_count = len(self.df)
        logging.debug(f"filter_by_columns end {name}, original_count: {original_count}, filtered_count: {filtered_count}")
        return self

    def select_by_name(self, name: str):
        logging.debug(f"select_by_name start {name}")
        original_count = len(self.df)
        self.df = self.df[self.df["FILTER"].str.contains(name, na=False)]
        filtered_count = len(self.df)
        logging.debug(f"select_by_name end {name}, original_count: {original_count}, filtered_count: {filtered_count}")
        return self

    def cut_filter(self, only_save_columns: list = None):
        logging.debug(f"cut_filter start")
        if 'FILTER' not in self.df.columns:
            logging.warning("FILTER not in columns")
            return self
        parsed_filter_list = []
        new_columns = set()
        for filter_entry in self.df["FILTER"]:
            entry_dict = {}
            for pair in filter_entry.split(";"):
                entry_dict[pair] = True
            new_columns.update(entry_dict.keys())
            parsed_filter_list.append(entry_dict)

        parsed_filter_df = pd.DataFrame(parsed_filter_list)
        parsed_filter_df.index = self.df.index
        self.df = pd.concat([self.df, parsed_filter_df], axis=1)

        if '.' in self.df.columns and '.' in new_columns:
            self.df.drop(columns=['.'], inplace=True)
        for col in self.df.columns:
            if only_save_columns is not None and col in new_columns and col not in only_save_columns:
                self.df.drop(columns=[col], inplace=True)

        logging.debug(f"new_columns: {new_columns}")
        logging.debug("cut_filter end ")
        return self

    def cut_info(self, only_save_columns: list = None):
        logging.debug(f"cut_info start")
        if 'INFO' not in self.df.columns:
            logging.warning("INFO not in columns")
            return self

        parsed_info_list = []
        new_columns = set()

        for info_entry in self.df["INFO"]:
            entry_dict = {}
            for pair in info_entry.split(";"):
                if "=" in pair:
                    key, value = pair.split("=", 1)
                    entry_dict[key] = value
                else:
                    # 處理布林標籤（無等號的鍵）
                    entry_dict[pair] = True
            new_columns.update(entry_dict.keys())
            parsed_info_list.append(entry_dict)

        parsed_info_df = pd.DataFrame(parsed_info_list)
        parsed_info_df.index = self.df.index
        self.df = pd.concat([self.df, parsed_info_df], axis=1)

        if '.' in self.df.columns and '.' in new_columns:
            self.df.drop(columns=['.'], inplace=True)
        for col in self.df.columns:
            if only_save_columns is not None and col in new_columns and col not in only_save_columns:
                self.df.drop(columns=[col], inplace=True)

        logging.debug(f"new_columns: {new_columns}")
        logging.debug("cut_info end ")
        return self

    def cut_format_sample(self, only_save_columns: list = None):
        logging.info("cut_format_sample start")
        if 'FORMAT' not in self.df.columns or 'SAMPLE' not in self.df.columns:
            logging.warning("FORMAT or SAMPLE not in columns")
            return self

        unique_formats = self.df['FORMAT'].unique()
        result_dfs = []
        new_columns = set()

        for fmt in unique_formats:
            subset = self.df[self.df['FORMAT'] == fmt]
            fmt_cols = fmt.split(':')
            new_columns.update(fmt_cols)
            subset_expanded = subset['SAMPLE'].str.split(':', expand=True)
            subset_expanded.columns = fmt_cols

            subset = pd.concat([subset, subset_expanded], axis=1)
            result_dfs.append(subset)

        self.df = pd.concat(result_dfs, axis=0).sort_index()
        if '.' in self.df.columns and '.' in new_columns:
            self.df.drop(columns=['.'], inplace=True)
        for col in self.df.columns:
            if only_save_columns is not None and col in new_columns and col not in only_save_columns:
                self.df.drop(columns=[col], inplace=True)

        logging.debug(f"new_columns: {new_columns}")
        logging.debug("cut_format_sample end")
        return self

    def check_variant_type(self):
        logging.debug("check_variant_type start")

        alt = self.df['ALT']
        ref = self.df['REF']

        # 打印 alt 或 ref 为 nan 的行
        nan_mask = alt.isna() | ref.isna()
        # print(nan_mask)
        if nan_mask.any():
            logging.error("发现 alt 或 ref 为 nan 的行：")
            logging.error(self.df[nan_mask])
        self.df = self.df[~nan_mask]
        alt = self.df['ALT']
        ref = self.df['REF']

        # 定義各種條件
        condition_multi_allelic = alt.str.contains(',')
        condition_insertion = ((alt.str.len() > 1) & (ref.str.len() == 1))
        condition_deletion = ((alt.str.len() == 1) & (ref.str.len() > 1))
        condition_snp = ((alt.str.len() == 1) & (ref.str.len() == 1))

        # 建立條件列表和對應的選擇值
        conditions = [
            condition_multi_allelic,
            condition_insertion,
            condition_deletion,
            condition_snp
        ]

        variant_type_choices = [
            'multi_allelic',
            'insertion',
            'deletion',
            'snp'
        ]

        variant_type_simple_choices = [
            'other',
            'indel',
            'indel',
            'snp'
        ]

        # 使用 numpy.select 來分配 'type_complex'
        self.df['type_complex'] = np.select(
            conditions,
            variant_type_choices,
            default='other'
        )

        # 使用 numpy.select 來分配 'type_simple'
        self.df['type_simple'] = np.select(
            conditions,
            variant_type_simple_choices,
            default='other'
        )

        self.df['type_complex'] = self.df['type_complex'].astype('category')
        self.df['type_simple'] = self.df['type_simple'].astype('category')

        logging.debug("check_variant_type end")
        return self

    def sort_chr_pos(self, dictionary: bool = False, inplace: bool = True):
        self.df = self.df.custom.sort_chr_pos(dictionary, inplace)
        return self


def find_range(range_df, pos_df, run_name):
    # logging.debug(f"{run_name} run")
    result = []
    if not range_df.empty:
        current_iterator = range_df.index[0]
        iterator_end = range_df.index[-1]
        pos_df["result"] = 0
        for index, row_pos in pos_df.iterrows():
            pos = row_pos["POS"]
            if not range_df.empty:
                for current in range(current_iterator, iterator_end + 1):
                    # https://samtools.github.io/hts-specs/VCFv4.2.pdf
                    # vcf POS 起始 1
                    # https://genome.ucsc.edu/FAQ/FAQformat.html#:~:text=chromStart%20%2D%20The%20starting,of%20a%20chromosome.
                    # bed chromStart 起始 0
                    # bed chromEnd 不包含結尾
                    range_start = range_df.loc[current, "POS1"] + 1
                    range_end = range_df.loc[current, "POS2"]
                    # debug use
                    # print(pos,range_start,range_end,current,current_iterator,iterator_end)
                    if pos > range_end:
                        continue
                    elif pos >= range_start and pos <= range_end:
                        pos_df.loc[index, "result"] = 1
                        break
                    elif pos < range_start:
                        break
                current_iterator = current
        result = pos_df.query("result == 1").index
    # logging.debug(f"{run_name} end ,result :{len(result)}")
    return result


class Diff_bed_vcf:
    def __init__(
        self,
        bed_df: pd.DataFrame = None,
        vcf_df: pd.DataFrame = None,
        parallel_multiple: int = 1,
    ):
        # logging.info(
        #     "Diff_bed_vcf---------------------------------------------------")
        self.bed_df = bed_df
        self.vcf_df = vcf_df
        if parallel_multiple < 1:
            parallel_multiple = 1
        self.parallel_multiple = parallel_multiple

    def find_range_muti_chrs(self, name: str = "INBED"):
        logging.info(f"find_range_muti_chrs start")
        pool = Pool(
            len(self.vcf_df["CHROM"].unique().tolist()) * self.parallel_multiple
        )
        result_list = []
        find_range_count = 0
        for chr, vcf_df2_chr in self.vcf_df.groupby("CHROM"):
            bed_df1_chr = self.bed_df.query("CHROM == @chr")[["POS1", "POS2"]]
            vcf_df2_chr_pos = vcf_df2_chr[["POS"]]
            num = len(vcf_df2_chr_pos) // self.parallel_multiple
            for i in range(self.parallel_multiple):
                if i != self.parallel_multiple - 1:
                    vcf_df2_chr_pos_cut_data = vcf_df2_chr_pos.iloc[
                        num * i: num * (i + 1)
                    ]
                else:
                    vcf_df2_chr_pos_cut_data = vcf_df2_chr_pos.iloc[num * i:]
                result_list.append(
                    pool.apply_async(
                        find_range,
                        args=(bed_df1_chr, vcf_df2_chr_pos_cut_data, f"{chr}-{i}")
                    )
                )
                find_range_count += len(vcf_df2_chr_pos_cut_data)
        pool.close()
        pool.join()
        self.vcf_df[name] = 0
        # self.vcf_df.loc[:, 'INBED'] = 0
        for result in result_list:
            result.wait()  # 應該是沒有甚麼用
            # print(result.get())
            self.vcf_df.loc[result.get(), name] = 1
        if find_range_count != len(self.vcf_df):
            logging.warning(f"find_range_muti_chrs end, find_range_count: {find_range_count}, vcf_df_count: {len(self.vcf_df)}")
        logging.debug(f"find_range_muti_chrs end")
        logging.debug(f"\n{self.vcf_df[name].value_counts()}")

        return self


def tag_germline(vcf_df, germline_df, use_alt=True):
    merge_cols = ["CHROM", "POS", "ALT"] if use_alt else ["CHROM", "POS"]
    vcf_df = vcf_df.merge(germline_df[merge_cols], on=merge_cols, how='left', indicator=True)

    name = "germline_tmp"
    if name not in vcf_df.columns and 'germline' in vcf_df.columns:
        vcf_df.rename(columns={'germline': name}, inplace=True)

        vcf_df['germline'] = (
            (vcf_df[name] == 1) |
            (vcf_df['_merge'] == 'both')
        ).astype(int)
        vcf_df.drop(columns=[name], inplace=True)
    else:
        vcf_df['germline'] = vcf_df['_merge'] == 'both'

    vcf_df.drop(columns=['_merge'], inplace=True)
    return vcf_df


def ans_tag(vcf_df, ans_df, use_alt=True):
    merge_cols = ["CHROM", "POS", "ALT"] if use_alt else ["CHROM", "POS"]
    vcf_df = pd.merge(vcf_df, ans_df[merge_cols], on=merge_cols, how="outer", indicator=True)

    mapping = {
        'both': 'tp',
        'left_only': 'fp',
        'right_only': 'fn'
    }
    vcf_df.custom.rename_merge_and_map(mapping, 'ans', inplace=True)
    return vcf_df


def get_germline():
    data_root_path = get_config_column("data_root_path")
    germline_1000g = BioinfoFileReader(file_path=f"{data_root_path}/PON/clairs-to_databases/1000g-pon.sites.vcf")
    germline_1000g.reader()
    germline_1000g.df = germline_1000g.df[['CHROM', 'POS', 'ALT']]

    germline_dbsnp = BioinfoFileReader(file_path=f"{data_root_path}/PON/clairs-to_databases/dbsnp.b138.non-somatic.sites.vcf")
    germline_dbsnp.reader()
    germline_dbsnp.df = germline_dbsnp.df[['CHROM', 'POS', 'ALT']]

    germline_gnomad = BioinfoFileReader(file_path=f"{data_root_path}/PON/clairs-to_databases/gnomad.r2.1.af-ge-0.001.sites.vcf")
    germline_gnomad.reader()
    germline_gnomad.df = germline_gnomad.df[['CHROM', 'POS', 'ALT']]

    germline_deepvariant = BioinfoFileReader(
        file_path=f"{data_root_path}/PON/clairs-to_databases/CoLoRSdb.GRCh38.v1.1.0.deepvariant.glnexus.af-ge-0.001.vcf")
    germline_deepvariant.reader()
    germline_deepvariant.df = germline_deepvariant.df[['CHROM', 'POS', 'ALT']]

    return germline_1000g.df, germline_dbsnp.df, germline_gnomad.df, germline_deepvariant.df


def tag_germline_all(vcf_df, germline_1000g, germline_dbsnp, germline_gnomad, germline_deepvariant, use_alt=[False, True, True, False]):
    vcf_df = tag_germline(vcf_df, germline_1000g, use_alt=use_alt[0])
    vcf_df = tag_germline(vcf_df, germline_dbsnp, use_alt=use_alt[1])
    vcf_df = tag_germline(vcf_df, germline_gnomad, use_alt=use_alt[2])
    vcf_df = tag_germline(vcf_df, germline_deepvariant, use_alt=use_alt[3])
    return vcf_df

def somatic_test(in_path):
    df = BioinfoFileReader(in_path, columns=["CHROM", "POS", "HIGHV", "LOWV", "DISAGREEV"]).reader().df
    # # df = df[df["HIGH"] == 0]
    # df["denominator"] = df["LOWV"] + df["DISAGREEV"]
    # df["numerator"] = df["LOWV"]
    # # df = df[df["denominator"] > 0]
    # # df = df[df["numerator"] > 0]
    # df["ratio_t"] = df["numerator"] / df["denominator"]
    # # df = df[df["ratio"] < 0.99]
    # # df["ratio"] = df["numerator"] / df["denominator"]


    return df

def purity_test(in_path, bed_path=None):
    df = BioinfoFileReader(in_path, columns=["CHROM", "POS", "SOMATIC", "HP1_REF", "HP1_ALT", "HP2_REF", "HP2_ALT", "CONFIDENCE"]).reader().df
    # df = df[df["CHROM"] != "chrX"]
    # df = df[df["CHROM"] != "chrY"]
    # df = df[df["CHROM"] == "chr2"]
    df = df[df["SOMATIC"] == 1]
    # df['VAF'] = (df['HP1_ALT'] + df['HP2_ALT']) / (df['HP1_REF'] + df['HP1_ALT'] + df['HP2_ALT'] + df['HP2_REF'])
    # print(len(df))
    # df = df[df["SOMATIC"] != 1]

    # run_root_path = get_config_column("run_root_path")
    # bed_path = "/big8_disk/data/HCC1395/SEQC2/High-Confidence_Regions_v1.2.bed"
    # # bed_path = f"{run_root_path}/test-cnv/one_loh.bed"
    # bed_df = BioinfoFileReader(bed_path).reader().df
    # if bed_df is not None:
    #     df = Diff_bed_vcf(bed_df, df, 3).find_range_muti_chrs("HIGH").vcf_df.query("HIGH == 1").copy()

    # df["check"]=df["CONFIDENCE"] > 0.75
    # df = df[(df["CONFIDENCE"] > 0.75) & (df["CONFIDENCE"] < 0.95)]
    # df = df[df["CONFIDENCE"] > 0.75]
    df["hp_set"] = df["HP1_ALT"] < df["HP2_ALT"]
    df["set1_alt"] = df.apply(
        lambda row: row["HP2_ALT"] if row["hp_set"] == 1
        else row["HP1_ALT"], axis=1
    )
    df["set2_alt"] = df.apply(
        lambda row: row["HP1_ALT"] if row["hp_set"] == 1
        else row["HP2_ALT"], axis=1
    )
    df["set1_ref"] = df.apply(
        lambda row: row["HP2_REF"] if row["hp_set"] == 1
        else row["HP1_REF"], axis=1
    )
    df["set2_ref"] = df.apply(
        lambda row: row["HP1_REF"] if row["hp_set"] == 1
        else row["HP2_REF"], axis=1
    )
    # df["het_ratio"] = (df["set1_alt"] + df["set1_ref"]) / (df["set1_alt"] + df["set1_ref"] + df["set2_ref"])
    # df = df[(df["het_ratio"] <= 0.55) & (df["het_ratio"] >= 0.45)]
    # df = df[(df["set1_alt"] + df["set1_ref"] >=5)]

    # df["denominator"] = df["set1_alt"] + df["set2_alt"] + df["set2_ref"]
    # df["numerator"] = df["set1_alt"] + df["set2_alt"]
    # df["denominator"] = df["set1_alt"] + df["set2_ref"]
    # df["numerator"] = df["set1_alt"]

    # df["max"] = df[["set2_ref", "set1_ref"]].apply(max, axis=1)
    # df["min"] = df[["set2_ref", "set1_ref"]].apply(min, axis=1)
    # df["numerator"] = df["max"] - df["min"] + df["set1_alt"] + df["set2_alt"]
    # df["denominator"] = df["set2_ref"] + df["set1_ref"] + df["set1_alt"] + df["set2_alt"]

    # df["test1"] = df["set2_ref"] + df["set1_alt"]
    # df["test2"] = df["set1_ref"] + df["set2_alt"]
    # df["denominator"] = df["test2"] + df["test1"]
    # df["numerator"] = df["test2"]

    df["denominator"] = df["set1_ref"] + df["set2_ref"]
    df["numerator"] = df[["set1_ref", "set2_ref"]].apply(max, axis=1)
    # df["numerator"] = df["set1_ref"]

    # df = df[df["denominator"] > 0]
    # df["ratio"] = df["numerator"] / df["denominator"]

    # df = df[df["HP2_REF"] > 1]
    # df["ratio"] = 1/df["HP2_REF"]
    # df["hp_ratio"] = (df["HP1_REF"] + df["HP1_ALT"]) / (df["HP1_REF"] + df["HP2_REF"] + df["HP1_ALT"] + df["HP2_ALT"])
    # df = df[(df["hp_ratio"] >= 0.40) & (df["hp_ratio"] <= 0.60)]

    df["all_ref"] = df["HP1_REF"] + df["HP2_REF"]
    df = df[df["all_ref"] > 1]
    df["ratio"] = df[["HP1_REF", "HP2_REF"]].apply(max, axis=1)/df["all_ref"]
    # df["ratio"] = df["HP1_REF"]/df["all_ref"]

    # df = df[df["ratio"] < 0.99]
    # print(df.query("ratio >0.98 "))
    # df = df[~(df["check"] & (df["ratio"] == 1))]
    return df


def phased_tag_ans(file_path, ans_vcf_path, ans_bed_path: str = None):
    file_vcf_path = file_path+".vcf"
    file_tmp_path = file_path+".body"
    file_tags_path = file_path+"_tags.vcf"
    calling_vcf_df = BioinfoFileReader(file_vcf_path).reader().df
    ans_vcf_df = BioinfoFileReader(ans_vcf_path).reader().df

    calling_somatic_df = calling_vcf_df.copy()
    calling_somatic_df["orig_index"] = calling_somatic_df.index
    calling_somatic_df = calling_somatic_df[calling_somatic_df['SAMPLE'].str.contains(r"0[|/]0", na=False)].copy()
    calling_somatic_df = pd.merge(calling_somatic_df, ans_vcf_df[["CHROM", "POS", "ALT"]], on=["CHROM", "POS", "ALT"], how="inner")
    if os.path.exists(ans_bed_path):
        ans_bed_df = BioinfoFileReader(ans_bed_path).reader().df
        calling_somatic_df = Diff_bed_vcf(ans_bed_df, calling_somatic_df, 3).find_range_muti_chrs("HIGH").vcf_df.query("HIGH == 1")
    # else:
    #     logging.warning(f"File {ans_bed_path} does not exist.")

    mask = calling_vcf_df.index.isin(calling_somatic_df["orig_index"])
    calling_vcf_df.loc[mask, "SAMPLE"] += ":1"
    calling_vcf_df.loc[mask, "FORMAT"] += ":ST"

    calling_vcf_df.to_csv(file_tmp_path, sep="\t", index=False, header=False)
    subprocess.run(f"grep '^#' {file_vcf_path} > {file_tags_path}", shell=True, check=True)
    subprocess.run(f"cat {file_tmp_path} >> {file_tags_path}", shell=True, check=True)
    return file_tags_path
