package mpi
import "fmt"
import "unsafe"

// struct ompi_communicator_t {};
// struct ompi_datatype_t {};
// #include <mpi.h>
import "C"

type Error struct {
	Code int;
}

type Status struct {
	Source, Tag, Error int;
}

func (err *Error) String() string {
	return fmt.Sprintf("MPI error with code %d", err.Code);
}

func Size() int {
	var size C.int;
	C.MPI_Comm_size(COMM_WORLD, &size);
	return int(size);
}

func Rank() int {
	var rank C.int;
	C.MPI_Comm_rank(COMM_WORLD, &rank);
	return int(rank);
}

func Init() int {
	return int(C.MPI_Init(nil, nil));
}

func Abort(errorCode int) {
	C.MPI_Abort(COMM_WORLD, C.int(errorCode));
}

func Finalize() {
	C.MPI_Finalize();
}

func SendInt32(data []int32, dest, tag int) int {
	// sizeof(MPI_INTEGER) == 4
	ret := C.MPI_Send(unsafe.Pointer(&data[0]), C.int(len(data)), MPI_INTEGER, C.int(dest), C.int(tag),
		COMM_WORLD);
	return int(ret);
}

func RecvInt32(data []int32, source, tag int) Status {
	var status C.MPI_Status;
	C.MPI_Recv(unsafe.Pointer(&data[0]), C.int(len(data)), MPI_INTEGER, C.int(source), C.int(tag),
		COMM_WORLD, &status);
	var s Status;
	s.Source = int(status.MPI_SOURCE);
	s.Tag = int(status.MPI_TAG);
	s.Error = int(status.MPI_ERROR);
	return s;
}

func SendFloat32(data []float32, dest, tag int) int {
        // sizeof(MPI_FLOAT) == 4
        ret := C.MPI_Send(unsafe.Pointer(&data[0]), C.int(len(data)), MPI_FLOAT, C.int(dest), C.int(tag),
                COMM_WORLD);
        return int(ret);
}

func RecvFloat32(data []float32, source, tag int) Status {
        var status C.MPI_Status;
        C.MPI_Recv(unsafe.Pointer(&data[0]), C.int(len(data)), MPI_FLOAT, C.int(source), C.int(tag),
                COMM_WORLD, &status);
        var s Status;
        s.Source = int(status.MPI_SOURCE);
        s.Tag = int(status.MPI_TAG);
        s.Error = int(status.MPI_ERROR);
        return s;
}


var COMM_WORLD = &C.ompi_mpi_comm_world;
var MPI_INTEGER = &C.ompi_mpi_integer;
var MPI_CHARACTER = &C.ompi_mpi_character;
var MPI_FLOAT = &C.ompi_mpi_float; // This is a float32 for us
var MPI_DOUBLE = &C.ompi_mpi_double; // float64

const (
 MPI_SUCCESS                   = 0; 
 MPI_ERR_BUFFER                = 1;
 MPI_ERR_COUNT                 = 2;
 MPI_ERR_TYPE                  = 3;
 MPI_ERR_TAG                   = 4;
 MPI_ERR_COMM                  = 5;
 MPI_ERR_RANK                  = 6;
 MPI_ERR_REQUEST               = 7;
 MPI_ERR_ROOT                  = 8;
 MPI_ERR_GROUP                 = 9;
 MPI_ERR_OP                    = 10;
 MPI_ERR_TOPOLOGY              = 11;
 MPI_ERR_DIMS                  = 12;
 MPI_ERR_ARG                   = 13;
 MPI_ERR_UNKNOWN               = 14;
 MPI_ERR_TRUNCATE              = 15;
 MPI_ERR_OTHER                 = 16;
 MPI_ERR_INTERN                = 17;
 MPI_ERR_IN_STATUS             = 18;
 MPI_ERR_PENDING               = 19;
 MPI_ERR_ACCESS                = 20;
 MPI_ERR_AMODE                 = 21;
 MPI_ERR_ASSERT                = 22;
 MPI_ERR_BAD_FILE              = 23;
 MPI_ERR_BASE                  = 24;
 MPI_ERR_CONVERSION            = 25;
 MPI_ERR_DISP                  = 26;
 MPI_ERR_DUP_DATAREP           = 27;
 MPI_ERR_FILE_EXISTS           = 28;
 MPI_ERR_FILE_IN_USE           = 29;
 MPI_ERR_FILE                  = 30;
 MPI_ERR_INFO_KEY              = 31;
 MPI_ERR_INFO_NOKEY            = 32;
 MPI_ERR_INFO_VALUE            = 33;
 MPI_ERR_INFO                  = 34;
 MPI_ERR_IO                    = 35;
 MPI_ERR_KEYVAL                = 36;
 MPI_ERR_LOCKTYPE              = 37;
 MPI_ERR_NAME                  = 38;
 MPI_ERR_NO_MEM                = 39;
 MPI_ERR_NOT_SAME              = 40;
 MPI_ERR_NO_SPACE              = 41;
 MPI_ERR_NO_SUCH_FILE          = 42;
 MPI_ERR_PORT                  = 43;
 MPI_ERR_QUOTA                 = 44;
 MPI_ERR_READ_ONLY             = 45;
 MPI_ERR_RMA_CONFLICT          = 46;
 MPI_ERR_RMA_SYNC              = 47;
 MPI_ERR_SERVICE               = 48;
 MPI_ERR_SIZE                  = 49;
 MPI_ERR_SPAWN                 = 50;
 MPI_ERR_UNSUPPORTED_DATAREP   = 51;
 MPI_ERR_UNSUPPORTED_OPERATION = 52;
 MPI_ERR_WIN                   = 53;
 MPI_ERR_LASTCODE              = 54;
 MPI_ERR_SYSRESOURCE           = -2;
)
