import MySQLdb as mdb
from sshtunnel import SSHTunnelForwarder


def query(sql_cmd):
    with SSHTunnelForwarder(
            ('airglowgroup.web.illinois.edu', 22),
            ssh_username='airglowgroup@web',
            ssh_private_key='/home/airglowgroup/.ssh/id_rsa',
            remote_bind_address=('127.0.0.1', 3306)
    ) as server:
        con = mdb.connect(host='127.0.0.1', db='airglowgroup_wp468',
                          port=server.local_bind_port, read_default_file="/home/airglowgroup/.my.cnf")
        cur = con.cursor()
        cur.execute(sql_cmd)
        rows = cur.fetchall()
        return rows


# site_id = site['sql_id']
# utc = pytz.utc  # Define timezones

# # Start and stop time of observations.
# if use_npz:
#     d0, dn = FPI_Results['sky_times'][0], FPI_Results['sky_times'][-1]
# else:
#     d = FPI.ReadIMG(sky_fns[0])
#     d0 = local.localize(d.info['LocalTime'])
#     d = FPI.ReadIMG(sky_fns[-1])
#     dn = local.localize(d.info['LocalTime'])

# startut = d0.astimezone(utc).strftime('%Y-%m-%d %H:%M:%S')
# stoput = dn.astimezone(utc).strftime('%Y-%m-%d %H:%M:%S')

# Open the database (see http://zetcode.com/databases/mysqlpythontutorial/ for help)
# Read the user and password from a file.
#        con = mdb.connect(host='webhost.engr.illinois.edu', db='airglowgroup_webdatabase', read_default_file="/home/airglow/.my.cnf")
#        cur = con.cursor()
# Send summary images to server
    # update the database
        # Send png. First find out if the entry is in there (i.e., we are just updating the png)

sql_cmd = "CREATE TABLE student (id INT AUTO_INCREMENT PRIMARY KEY, name VARCHAR(255));"
#sql_cmd_2 = "INSERT INTO student (name) VALUES(%s)", ("Cristina",)
# sql_cmd = 'SELECT id FROM DataSet WHERE Site = %d and Instrument = %d and StartUTTime = \"%s\"' % (site_id, db_id, startut)
rows = query(sql_cmd)
#            cur.execute(sql_cmd)
#            rows = cur.fetchall()
#rr = query(sql_cmd_2)
#     log_fn = db_log_stub + instrsitedate + '_log.log'
#     if len(rows) == 0: # Create the entry
#         sql_cmd = 'INSERT INTO DataSet (Site, Instrument, StartUTTime, StopUTTime, SummaryImage, LogFile) VALUES(%d, %d, \"%s\", \"%s\", \"%s\", \"%s\")' % (site_id, db_i>
#         logfile.write(sql_cmd + '\n')
#         query(sql_cmd)
# #                cur.execute(sql_cmd)
#     else: # Entry exists. Update it.
#         sql_cmd = 'UPDATE DataSet SET SummaryImage=\"%s\",LogFile=\"%s\" WHERE id = %d' % (db_image_stub + fn, log_fn, rows[0][0])
#         # sql_cmd = 'UPDATE DataSet SET SummaryImage=\"%s\",LogFile=\"%s\" WHERE Site = %d and Instrument = %d and StartUTTime = \"%s\"' % (db_image_stub + fn, log_fn, si>
# #                cur.execute(sql_cmd)
#         logfile.write(sql_cmd + '\n')
#         query(sql_cmd)
