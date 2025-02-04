import pm4py
import pandas as pd

def import_log(address):
    case_id_name = 'case:concept:name'
    timestamp_name = 'time:timestamp'
    acyivity_name = 'concept:name'

    log = pm4py.read_xes(str(address), variant="rustxes")
    case_table, event_table, map_info = log_to_tables(log, parameters={'case_id': case_id_name,
                                                                                        'timestamp': timestamp_name,
                                                                                        'activity_name': acyivity_name})
    case_table.to_csv('output_files/out.csv', index=False)
    event_table.to_csv('output_files/out_event.csv', index=False)


    return len(case_table),case_table.select_dtypes(include=['number']).columns


def log_to_tables(df, parameters):

    # Extract parameters
    case_id_name = parameters['case_id']
    timestamp_name = parameters['timestamp']
    activity_name = parameters['activity_name']

    # Constants
    time_unit = 24 * 3600  # Time unit in seconds (1 day)

    # Output column names
    output_case_id_name = 'case_id'
    output_timestamp_name = 'timestamp'
    output_activity_name = 'activity_name'

    # Identify case attributes (columns with unique values per case)
    dfunique = df.groupby(case_id_name).nunique()
    case_attributes = set(dfunique.columns[dfunique.max() <= 1])

    # Build the case table
    grouped = df.groupby(case_id_name)
    case_table = grouped.first()[list(case_attributes)]
    case_table['trace'] = grouped[activity_name].agg(lambda x: tuple(x.astype(str)))
    case_table['duration'] = (grouped[timestamp_name].max() - grouped[
        timestamp_name].min()).dt.total_seconds() // time_unit
    case_table['n_events'] = grouped.size()
    case_table = case_table.reset_index().rename(columns={case_id_name: output_case_id_name})

    # Build the event table
    non_attributes = {'row_num'}
    event_attributes = set(df.columns) - case_attributes - non_attributes
    event_table = df[list(event_attributes)]
    event_table[timestamp_name] = pd.to_datetime(event_table[timestamp_name], utc=True)

    # Sort events by timestamp (and EventOrder if available)
    sort_columns = [timestamp_name]
    if 'EventOrder' in event_table.columns:
        sort_columns.append('EventOrder')
    event_table = event_table.sort_values(sort_columns)

    # Rename columns for the event table
    event_table = event_table.rename(columns={
        case_id_name: output_case_id_name,
        timestamp_name: output_timestamp_name,
        activity_name: output_activity_name
    })

    # Mapping for standardized column names
    mapping = {
        'case_id': output_case_id_name,
        'timestamp': output_timestamp_name,
        'activity_name': output_activity_name
    }

    return case_table, event_table, mapping