module mod_sys

use iso_c_binding

interface
    function rename_wrapper(filein, fileout) bind(C, name="rename")
    import :: c_char, c_int
    integer(c_int) :: rename_wrapper
    character(kind=c_char) :: filein(*), fileout(*)
    endfunction rename_wrapper
endinterface

endmodule mod_sys
