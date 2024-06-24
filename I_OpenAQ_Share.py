import pandas as pd
import os
import shutil
import glob

from datetime import datetime, timedelta

# find the working direc
def find_root_dir(debug_mode):
    currentdir = os.getcwd()
     # For mac OS
    if 'Volumes' in currentdir and debug_mode == 0:
        direc='/Volumes/rvmartin/Active/SPARTAN-shared/'
        
    elif  'Volumes' in currentdir and debug_mode == 1:
        direc='/Volumes/rvmartin/Active/haihuizhu/6.SPARTAN/' 
    
    elif 'storage1/fs1' in currentdir and debug_mode == 0: 
        direc='/storage1/fs1/rvmartin/Active/SPARTAN-shared' 

    elif 'storage1/fs1' in currentdir and debug_mode == 1:
        direc='/storage1/fs1/rvmartin/Active/haihuizhu/6.SPARTAN/' 
        
    # for Windows OS
    elif '\\storage1.ris.wustl.edu'  in currentdir:
        direc= '\\\\storage1.ris.wustl.edu\\rvmartin\\Active\\SPARTAN-shared\\'
    elif 'Z:' in currentdir:
        direc= '\\\\storage1.ris.wustl.edu\\rvmartin\\Active\\SPARTAN-shared\\'
    else:
        # raise ValueError('Cannot identify the root directory') 
        direc='/Volumes/rvmartin/Active/SPARTAN-shared/'# test only
    return direc
    

# Define the file path
debug_mode = 0;
direc = find_root_dir(debug_mode);
direc_in = f"{direc}/Public_Data/Time-resolved_PM2.5/Data_sharing/";
direc_out = f"{direc}/ExternalShare/ForOpenAQ/";

# read site sampling data
site_info = pd.read_excel(f"{direc}/Site_Sampling/Site_details.xlsx")

# Load OpenAQ recieved data
localtimerecord = f"{direc}/Public_Data/Scripts/UtilityFunctions/Received_by_OpenAQ_localtime.csv"
if os.path.exists(localtimerecord):
    openaq_rec = pd.read_csv(localtimerecord, parse_dates=['datetime_first','datetime_last']) 

else:
    openaq_rec = pd.read_csv(f"{direc}/Public_Data/Scripts/UtilityFunctions/Received_by_OpenAQ.csv", parse_dates=['datetime_first','datetime_last'])  
    # Create a dictionary to map city names to their time zones
    city_to_timezone = dict(zip(site_info['City'], site_info['Time_zone_GMT']))
    Ins_to_timezone = dict(zip(site_info['Host_Institute'], site_info['Time_zone_GMT']))
    # Convert datetimes in received_data to the corresponding time zones
    for index, row in openaq_rec.iterrows():
        city_name = row['name']
        if city_name in city_to_timezone:
            timezone = city_to_timezone[city_name] 
        elif city_name in Ins_to_timezone:
            timezone = Ins_to_timezone[city_name] 
        openaq_rec.loc[openaq_rec['name'].str.contains(city_name), 'datetime_first'] += timedelta(hours=timezone)
        openaq_rec.loc[openaq_rec['name'].str.contains(city_name), 'datetime_last'] += timedelta(hours=timezone)

        if city_name in Ins_to_timezone:
            new_name = site_info.loc[site_info['Host_Institute'].str.contains(city_name), 'City']
            openaq_rec.loc[openaq_rec['name'].str.contains(city_name), 'name'] = new_name.iloc[0] 
        else:
            print(f'{city_name} not found!')
    # save this file for future use
    openaq_rec.to_csv(localtimerecord, index=False)


# delete files in direc_out (Open_AQ doesn't want repeated data) 
csv_files = glob.glob(f'{direc_out}/*.csv')
print(f"Deleting data in {direc_out}")
for file_path in csv_files:
    try:
        os.remove(file_path)
    except Exception as e:
        print(f"Error deleting file {file_path}: {e}")


# loop througth site we have to get the time res est PM2.5
for index, row in site_info.iterrows():
    site_code = row['Site_Code']
    site_city = row['City']
    site_city_short = row['City'].replace(' ', '')
    timeres_in = f'{direc_in}{site_code}_{site_city_short}_DailyEstimatedPM25.csv'
    if os.path.exists(timeres_in):
        timeres_out = f'{direc_out}{site_code}_{site_city_short}_DailyEstimatedPM25.csv'
        timeres = pd.read_csv(timeres_in, skiprows=1)
        
        # if site in openaq dataset, find the right timeframe, save selected data to direct_out
        if openaq_rec['name'].str.contains(site_city).any():
            # locate the datetime_last in timeres
            timeres['combined_datetime']= pd.to_datetime(timeres['Year_local'].astype(str) + '-' +
                                                      timeres['Month_local'].astype(str) + '-' +
                                                      timeres['Day_local'].astype(str)
                                                      )  
            date_last = openaq_rec.loc[openaq_rec['name'].str.contains(site_city), 'datetime_last'].iloc[0]
            filtered_timeres = timeres[timeres['combined_datetime'] > date_last.tz_convert(None)]
            filtered_timeres = filtered_timeres.drop('combined_datetime', axis=1)
            # copy data after to a new file in direct_out
            filtered_timeres.to_csv(timeres_out, index=False)
            
        # if site not in openaq dataset, add the entire file to direct_out and add info to openaq_rec
        else:    
            # copy file
            shutil.copy(timeres_in, timeres_out)  
            # add to openaq_rec
            comb_time = pd.to_datetime(timeres['Year_local'].astype(str) + '-' +
                                       timeres['Month_local'].astype(str) + '-' +
                                       timeres['Day_local'].astype(str))
            new_row = {
                'name': site_city,
                'datetime_first': comb_time.iloc[0],
                'datetime_last':comb_time.iloc[-1]
            }
            new_row_df = pd.DataFrame([new_row])
            df = pd.concat([openaq_rec, new_row_df], ignore_index=True)   
            # Save the openaq_rec for future use
            openaq_rec.to_csv(localtimerecord, index=False)
            
        # change file permission # doesn't seem to work
        os.chmod(timeres_out, 0o644)
        
print('Done')