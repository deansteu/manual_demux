import pandas as pd
import boto3

class DemultiplexInstance:
    def __init__(self, sample_sheet, config:dict):
        config = {
            'output_dir': 'path/to/output'
            ''
        
        }

        self.config = config
        self.s3 = boto3.client('s3')
        self.demux_class = 'bulk or sc'
        
    def run_demux(self):
        pass

    def upload_project_fastqs(self):
        bucket = None
        if self.demux_class == 'bulk':
            bucket = 'cumc-ngs-u1-fastq/NovaSeq'
        elif self.demux_class == 'sc':
            bucket = 'cumc-scac-u1-fastq/'

    def upload_undetermined_fastqs(self):
        bucket = 'cumc-scac-u1-undetermined'
        key = '<RUN_ID>'

    def upload_lane_barcode_html(self):
        bucket = 'cumc-lane-barcode'




class Run:
    def __init__(self):
        pass

    def run_demux_instances():
        pass

    def scc_sample_sheets(self, df: pd.DataFrame, header_map: dict):
        scc_output_file = 'scc.csv'
        scc_df = df[[header_map['Lane'], header_map['Sample'], header_map['Index']]].copy()
        scc_df.columns = ['Lane', 'Sample', 'Index']
        scc_df.to_csv(scc_output_file, index=False)
        print(f'Data written to: {scc_output_file}')



    def reverse_complement(self, sequence: str) -> str:
        seq_rev = sequence[::-1]
        seq_rev_comp = ""
        for char in seq_rev:
            if char == 'A':
                seq_rev_comp += 'T'
            elif char == 'T':
                seq_rev_comp += 'A'
            elif char == 'G':
                seq_rev_comp+= 'C'
            elif char == 'C':
                seq_rev_comp += 'G'
            else:
                raise KeyError
        
        return seq_rev_comp

    def write_sample_sheet(self, df: pd.DataFrame, index: str, header_map: dict):
        data = {
            'Lane': df[header_map['Lane']],
            'Sample_ID': df[header_map['Sample']],
            'Sample_Name': df[header_map['Sample']],
            'Sample_Plate': '',
            'Sample_Well': '',
            'I7_Index_ID': df[header_map['Index']],
            'index': df['index'],  # Use the 'i7' column from the index DataFrame
            'I5_Index_ID': df[header_map['Index']],
            'index2': df['index2'].apply(lambda x: self.reverse_complement(x)),  # Use the 'i5' column from the index DataFrame
            'Sample_Project': df[header_map['Project']], #df['Project ID'],
            'Description': ''
        }
        output_df = pd.DataFrame(data)

        output_file = 'index_'+index+'.csv'
        with open(output_file, 'w', newline='') as f:
            # Write the header
            f.write('[Header]\n')
            f.write('EMFileVersion,4\n\n')
                        
            # Write the reads section
            f.write('[Reads]\n')
            f.write('150\n')
            f.write('150\n\n')
                        
            # Write the data section
            f.write('[Data]\n')
                        
            # Write the DataFrame without extra newlines
            output_df.to_csv(f, index=False)

    def bulk_sample_sheets(self, df: pd.DataFrame, header_map: dict):
        # index 16 or index 20
        # should clean this up, shouldn't have to do this, fix in the merge
        df['INDEX'] = df['INDEX_x']
        df.drop(['INDEX_x', 'INDEX_y'], axis=1, inplace=True)
        print(len(df))
        df[['index', 'index2']] = df['INDEX'].str.split('-', expand=True)

        df_20 = df[df['INDEX'].str.len() == 21]
        df_16 = df[df['INDEX'].str.len() == 17]
        # other index lengths?


        print(df_20)
        print(df_16)
        
        if not df_16.empty:
            self.write_sample_sheet(df_16, '16', header_map)

        if not df_20.empty:
            self.write_sample_sheet(df_16, '20', header_map)


    def read_ingest_sheet(self, ingest_sheet: str):
        """The original sheet pulled from HIQUEUE or V2_SINGLE_CELL
        """
        # read the full ingest sheet (would include multiple types of index lengths)
        # we will separate out whole ingest sheet into multiple sample sheets based on the index length and protocols
        if ingest_sheet.endswith('.tsv'):
            ingest_df = pd.read_csv(ingest_sheet, sep='\t', header=0)
        elif ingest_sheet.endswith('.csv'):
            ingest_df = pd.read_csv(ingest_sheet, header=0)
        else:
            print('unrecognized file extenstion in {}. Accepted formats are .csv and .tsv'.format(ingest_sheet))
            exit(1)
        

        dual_index_df = pd.read_csv('./index_sequences/Dual_Index_Kit_TT_TS_NN_NT_TN.csv')
        single_index_df = pd.read_csv('./index_sequences/Single_Index_NA.csv')
        bulk_index_df = pd.read_csv('./index_sequences/NexteraXT.csv')


        scc_headers = ['Project ID', 'Sample', 'Index', 'lane', 'Protocol']
        ngs_bulk_headers = ['LANE', 'JOB ID', 'SAMPLE ID', 'INDEX ID']
        demux_type = None
        if pd.Series(scc_headers).isin(ingest_df.columns).all():
            header_map = {
                'Lane': 'lane',
                'Index': 'Index',
                'Project': 'Project ID',
                'Sample': 'Sample'
            }
            demux_type='scc'
        elif pd.Series(ngs_bulk_headers).isin(ingest_df.columns).all():
            header_map = {
                'Lane': 'LANE',
                'Index': 'INDEX ID',
                'Project': 'JOB ID',
                'Sample': 'SAMPLE ID'
            }
            demux_type='bulk' 
        else:
            raise pd.errors.DataError(f'Sample Sheet {ingest_sheet} must have all column headers for one of the following\nFor Single Cell Data: {scc_headers}\nFor Bulk Data: {ngs_bulk_headers}\nHeaders were: {list(ingest_df.columns)}')

        print(ingest_sheet)
        print(ingest_df)
        print(ingest_df['INDEX ID'])
        # split the ingest sheet into the index length 
        # options
        # split index 8
        print(header_map['Index'])
        print(len(ingest_df['INDEX ID'].unique()))
        scc_df = pd.merge(ingest_df, dual_index_df, left_on=header_map['Index'], right_on='index_name', how='inner')
        scc_atac_df = pd.merge(ingest_df, single_index_df, left_on=header_map['Index'], right_on='index_name', how='inner')
        bulk_df = pd.merge(ingest_df, bulk_index_df, left_on=header_map['Index'], right_on='ID', how='inner')
        # print('scc df')
        # print(scc_df.columns)
        # print(scc_df)

        # print('scc atac df')
        # print(scc_atac_df.columns)
        # print(scc_atac_df)

        # print('bulk df')
        # print(bulk_df.columns)
        # print(bulk_df)

        results = pd.concat([scc_df, scc_atac_df, bulk_df])
        print('recombined df')
        print(results)
        print(results['LANE'].unique())
        print(len(results['INDEX ID'].unique()))
        missing_rows = ingest_df[~ingest_df[header_map['Index']].isin(results[header_map['Index']])]

        # print(missing_rows)
        non_indexed_df = missing_rows.copy()

        # print(f"original len: {len(ingest_df)}\ncombined len: {len(results)}:\nmissing: {len(non_indexed_df)}")

        # can we assume that if they don't appear in dual-index, single-index, or nextera, 
        # they are bulk (or should we treat them as something "separate")?
        if not scc_df.empty:
            self.scc_sample_sheets(scc_df, header_map)
        
        if not scc_atac_df.empty:
            self.scc_sample_sheets(scc_atac_df, header_map)

        if not bulk_index_df.empty:
            self.bulk_sample_sheets(bulk_df, header_map)
        
        if not non_indexed_df.empty:
            print("MUST DO SOMETHING WITH NON_INDEXED_DF!!!")

        # split index 8
        # split index 16
        # split index 20
        
        # single cell

        # bulk



    def create_demux_sample_sheet():
        pass


if __name__ == '__main__':
    print('hello')
    demux = Run()
    demux.read_ingest_sheet('./hq_ingest.tsv')