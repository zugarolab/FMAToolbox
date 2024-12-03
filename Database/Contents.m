% Database and batch processing functions for FMAToolbox.
%
% These functions can be used to perform batch data processing and store data
% and figures into one or more databases.
%
% Batch Processing
%
% Batch processing is useful if you need to repeatedly perform a given analysis
% on different data sets. This can be achieved using <a href="StartBatch">StartBatch</a>, which will
% repeatedly run a given 'batch' function on each data set. <a href="StartBatch">StartBatch</a> provides
% a number of advanced features such as error resilience, error logging, delayed
% processing, etc. The 'batch' function typically performs computations and
% stores the results in a database using <a href="DBAddVariable">DBAddVariable</a> and/or <a href="DBAddFigure">DBAddFigure</a>. Once
% the batch is completed, the results can be retrieved from the database using
% <a href="DBGetValues">DBGetValues</a>. Different sets of results can be combined using <a href="DBMatchValues">DBMatchValues</a>.
%
% Database Management
%
% Databases can be private and run on the local machine, or central databases
% running on a remote server and shared between many different users. Data and
% figures can be tagged (experiment ID, name, comments) so you can later search
% specific data subsets. Programs and parameters used to generate the data can
% also be stored for reference. Figures can be stored as matlab fig files
% and/or PNG images, and be later exported as HTML galleries which are easy to
% share with other investigators.
%
% Databases are optimized to store huge numbers of small items. However, Matlab
% variables or figures can be fairly large (think of a spectrogram recorded over
% serveral hours), and cannot be efficiently handled by a database. The general
% solution is to save such data to the file system, and keep a reference to their
% storage location in the database. In FMAToolbox, for the sake of consistency, all
% data (small or large) are thus saved to the file system. This will be transparently
% handled by the database functions. You must simply indicate the directory where
% Matlab files should be stored, and make sure it is shared and writable by all
% users who are expected to use the database. This can typically be done by adding
% a line in the startup.m file, for instance:
%
%   global SETTINGS;
%   SETTINGS.db.storage.externalPath = '/mnt/matlab-databases';
%
% An option is provided for specific use cases to store all or part of the data
% directly into the database. This is achieved by setting a non-zero upper limit
% to the size of variables and figures that will be stored in the database. Beyond
% this limit, the data will be stored as described above. To customize this limit,
% add one line in the startup.m file, for instance:
%
%   global SETTINGS;
%   SETTINGS.db.storage.maxInternalSize = 10*1024*1024;  % limit of 10MiB
%
% Databases use a number of fields to describe data (see e.g. <a href="matlab: help DBAddVariable">DBAddVariable</a>). The
% maximum lengths can be customized in the startup.m file:
%
%   global SETTINGS;
%   SETTINGS.db.fields.eid = 100; % string of max 100 characters
%   SETTINGS.db.fields.name = 100;
%   SETTINGS.db.fields.comments = 300;
%   SETTINGS.db.fields.parameters = 100;
%
% Obviously, changing these default values will not affect databases that have
% already been created.
%
% Notes
%
%  The database backend is MariaDB/MySQL and the interface is provided by the <a href="http://sourceforge.net/projects/mym/">mYm</a>
%  toolbox by Y. Maret (EPFL, Lausanne, Switzerland).
%
%  Although not required, installing <a href="http://sites.google.com/site/oliverwoodford/software/export_fig">export_fig</a> by O. Woodford (Oxford
%  University, England) will yield better PNG images.
%
% Batch processing
%
%   StartBatch            - Start a new batch job.
%   ShowBatch             - Show data sets in a batch job.
%   GetBatch              - Get batch job output.
%   BatchInfo             - Get batch job information.
%   CancelBatch           - Cancel batch job.
%   CleanBatches          - Delete completed batch jobs from memory.
%   DebugBatch            - Assign variables to help debug a batch job.
%   Store                 - Store variable in memory to avoid redundant computations.
%   Recall                - Recall variable from memory to avoid redundant computations.
%
% Database management.
%
%   DBCreate              - Create a new database (requires admin rights)
%   DBCreateTables        - Create tables for variables and figures
%   DBConnect             - Connect to the database server
%   DBExternalStoragePath - Get path for database external storage space
%   DBUse                 - Set (or determine) current database
%   DBList                - List existing databases
%   DBDuplicate           - Duplicate database.
%   DBMerge               - Merge databases.
%   DBRemove              - Remove database.
%
% Storing data/figures.
%
%   DBAddFigure           - Store a figure (+ tags, programs, etc.)
%   DBAddVariable         - Store a variable (+ tags, programs, etc.)
%   DBRemoveFigures       - Remove all figures that match given criteria.
%   DBRemoveVariables     - Remove all variables that match given criteria.
%
% Searching data/figures.
%
%   DBGetFigures          - Get all figures that match given criteria
%   DBGetValues           - Group values for all variables that match given criteria
%   DBMatchValues         - Match values for two different sets of variables
%   DBGetVariables        - Get all variables that match given criteria
%   DBListFigures         - List figures matching certain criteria.
%   DBListVariables       - List variables matching certain criteria.
%
% Exporting figures.
%
%   DBExportGallery       - Create an HTML gallery from a set of figures
%
