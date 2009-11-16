package main

import ( 
    "container/vector"; 
    "flag";
    "fmt";
    "go/ast";
    "go/parser";
    "go/token";
    "os";
    "path";
    "strings"; )

var sourceExt = ".go";
var objExt = ".?";
var depExt = ".d";

type importCollector interface {
    Push( string );
}


func collectImports( filename string, collector importCollector ) bool {
    parsedFile, err := parser.ParseFile( filename, nil, parser.ImportsOnly );
    if err != nil {
        fmt.Printf("Failed to parse '%v' %s.\n", filename, err);
        return false
    }
    
    for _, decl := range parsedFile.Decls {
        genDecl, ok := decl.(*ast.GenDecl);
        if !ok { continue; }
        if genDecl.Tok != token.IMPORT { continue; }
        for _, spec := range genDecl.Specs {
            importSpec, ok := spec.(*ast.ImportSpec);
            if !ok { continue; }
            importPath := "";
            for _, pathItem := range importSpec.Path {
                importPath += string(pathItem.Value);
            }
            collector.Push( importPath );
        }
    }
    
    return true;
}


func writeDepFile(sourcefileName string) bool {
    if path.Ext( sourcefileName ) != sourceExt {
        fmt.Printf("Error %v is missing the sourcefile suffix %v.\n", sourcefileName, sourceExt );
        fmt.Printf("It will be skipped.\n");
        return false;
    }

    importNames := vector.NewStringVector(0);
    if !collectImports( sourcefileName, importNames ) {
        return false;
    }

    dirname, basename := path.Split( sourcefileName );
    basename = basename[0:len(basename)-len( sourceExt)];
    depfileName := basename + depExt;    
        
    depFile, err := os.Open( path.Join( dirname, depfileName), os.O_WRONLY | os.O_CREATE | os.O_TRUNC, 0666 );
    if err != nil {
        fmt.Printf("Failed to open output file %v.\n", depfileName );
        return false;
    }
    depFile.Write( strings.Bytes(basename + objExt + " " + depfileName + " : " + basename + sourceExt ) );
    for _, importName := range importNames.Data() {
        if strings.HasPrefix( importName, "\"./" ) && len(importName) >= 5 {
            relativeName := importName[3:len(importName)-1];
            depFile.Write( strings.Bytes( " " + relativeName + objExt ) );
        }
    }
    depFile.Write( strings.Bytes( "\n" ) );
    depFile.Close();
    
    return true;
}


func main() {
    switch GOARCH := os.Getenv("GOARCH"); GOARCH {
        case "amd64": objExt = ".6";
        case "386": objExt = ".8";
        case "arm": objExt = ".5";
    }
    
    flag.StringVar( &sourceExt, "sourceExt", sourceExt, "The file extension for source files." );
    flag.StringVar( &objExt, "objExt", objExt, "The file extension for object files." );
    flag.StringVar( &depExt, "depExt", depExt, "The file extension for generated dependency files." );
    flag.Parse();
    
    if flag.NArg() == 0 {
        fmt.Printf("usage: depscan [-sourceExt <srcExtension>] [-objExt <objExtension>] [-depExt <depExtension>] sourcefile ...\n");
        os.Exit(-1);
    }
    
    ok := true;
    for i := 0; i < flag.NArg(); i++ {
        ok = ok && writeDepFile( flag.Arg(i) );     
    }
    if ok { os.Exit(0); }
    else { os.Exit(-1); }
}
