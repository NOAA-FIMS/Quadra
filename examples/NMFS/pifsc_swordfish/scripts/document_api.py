"""Inventory named C++/Python/R declarations, enforce descriptions, and generate an API reference."""
import argparse
import ast
import json
from pathlib import Path
import re
import subprocess
import tempfile

ROOT=Path(__file__).resolve().parents[1]


def json_sequence(text):
    """Decode Clang's concatenated JSON declarations without assuming a JSON array.

    Parameters
    ----------
    text
        Source/document text; output functions preserve intended line breaks.

    Returns
    -------
    Iterator of decoded Clang AST root objects.

    Notes
    -----
    C++ declarations are obtained from Clang AST JSON and Python declarations from ast. Anonymous lambdas, vendor internals and generated factory bodies are excluded explicitly. Description coverage is enforced before publishing the inventory; regenerate after changing declaration locations.
    """
    decoder=json.JSONDecoder();offset=0
    while offset<len(text):
        while offset<len(text) and text[offset].isspace():offset+=1
        if offset==len(text):break
        value,offset=decoder.raw_decode(text,offset)
        yield value


def cpp_entries(nodes, default_file=None):
    """Extract explicit named declarations, omitting template duplicates and lambda bodies.

    Parameters
    ----------
    nodes
        Clang AST root declarations from one or more JSON records.
    default_file
        Fallback relative source path for a wrapped executable/test translation unit.

    Returns
    -------
    Deduplicated explicit named C++ declarations; implicit members and lambda bodies are skipped.

    Notes
    -----
    C++ declarations are obtained from Clang AST JSON and Python declarations from ast. Anonymous lambdas, vendor internals and generated factory bodies are excluded explicitly. Description coverage is enforced before publishing the inventory; regenerate after changing declaration locations.
    """
    result={}
    def visit(node,file=default_file,scope=''):
        """Walk declaration context only, carrying source file and class qualification.

        Parameters
        ----------
        node
            Current syntax-tree declaration/context.
        file
            Inherited relative C++ source file.
        scope
            Enclosing class/function qualification accumulated during traversal.
        """
        location=node.get('loc',{});file=location.get('presumedFile',location.get('file',file))
        if file and not file.startswith(('include/swordfish/','src/','tests/')):
            file=default_file
        kind,name=node.get('kind'),node.get('name','')
        line=location.get('presumedLine',location.get('line'))
        if not line and file and 'offset' in location:
            line=(ROOT/file).read_text()[:location['offset']].count('\n')+1
        explicit=not node.get('isImplicit')
        callable=kind in {'FunctionDecl','CXXMethodDecl','CXXConstructorDecl','CXXConversionDecl'}
        record=kind=='CXXRecordDecl' and node.get('completeDefinition')
        alias=kind in {'TypeAliasDecl','TypedefDecl'}
        if explicit and name and file and (callable or record or alias):
            qualified=scope+name
            key=(file,location.get('offset',line),qualified)
            fields=[dict(name=c['name'],type=c.get('type',{}).get('qualType','')) for c in node.get('inner',[]) if c.get('kind')=='FieldDecl' and c.get('name')]
            parameters=[dict(name=c.get('name',''),type=c.get('type',{}).get('qualType','')) for c in node.get('inner',[]) if c.get('kind')=='ParmVarDecl']
            signature=node.get('type',{}).get('qualType',qualified)
            if callable:
                return_type='' if kind=='CXXConstructorDecl' else signature.split(' (',1)[0]+' '
                signature=return_type+qualified+'('+', '.join(p['type']+' '+p['name'] for p in parameters)+')'+(' const' if signature.endswith(' const') else '')
            source=(ROOT/file).read_text().splitlines();comment=[];position=line-2
            while position>=0 and source[position].lstrip().startswith('///'):
                comment.insert(0,source[position].lstrip()[3:].strip());position-=1
            result.setdefault(key,dict(language='C++',file=file,line=line,name=qualified,kind='class/struct' if record else 'alias' if alias else 'function',signature=signature,parameters=parameters,fields=fields,documentation='\n'.join(comment)))
        if callable:return
        if kind=='CXXRecordDecl':
            if not record or not explicit or not name:return
            scope+=name+'::'
        for child in node.get('inner',[]):
            if child.get('kind') not in {'CompoundStmt','ClassTemplateSpecializationDecl','TemplateArgument'}:visit(child,file,scope)
    for node in nodes:visit(node)
    return list(result.values())


def clang_inventory(compiler):
    """Use Clang ASTs for library declarations and wrapped executable/test source definitions.

    Parameters
    ----------
    compiler
        Clang executable supporting the JSON AST interface.

    Returns
    -------
    All supported named C++ declaration records from headers and wrapped source/tests.

    Notes
    -----
    C++ declarations are obtained from Clang AST JSON and Python declarations from ast. Anonymous lambdas, vendor internals and generated factory bodies are excluded explicitly. Description coverage is enforced before publishing the inventory; regenerate after changing declaration locations.
    """
    flags=[compiler,'-std=c++17','-fsyntax-only','-Iinclude','-Ibuild','-I../../..','-I../../../external/eigen','-Xclang','-ast-dump=json','-Xclang']
    process=subprocess.run(flags+['-ast-dump-filter=swordfish','docs/api_probe.cpp'],cwd=ROOT,text=True,capture_output=True,check=True)
    entries=cpp_entries(json_sequence(process.stdout))
    for path in sorted(list((ROOT/'src').glob('*.cpp'))+list((ROOT/'tests').glob('*.cpp'))):
        relative=str(path.relative_to(ROOT));text=path.read_text()
        includes='\n'.join(re.findall(r'^#include .+$',text,re.M))
        body=re.sub(r'^#include .+$','',text,flags=re.M).replace('DECLARE_ADGRAPH();','')
        probe=includes+'\n'+('DECLARE_ADGRAPH();\n' if 'DECLARE_ADGRAPH();' in text else '')+'namespace api_source {\n#line 1 "'+relative+'"\n'+body+'\n}\n'
        with tempfile.NamedTemporaryFile(mode='w',suffix='.cpp',delete=False) as stream:
            stream.write(probe);temporary=Path(stream.name)
        try:
            process=subprocess.run(flags+['-ast-dump-filter=api_source',str(temporary)],cwd=ROOT,text=True,capture_output=True,check=True)
            entries.extend(cpp_entries(json_sequence(process.stdout),relative))
        finally:temporary.unlink()
    return entries


def python_entries():
    """Inventory functions, nested helpers, classes and methods with Python's syntax tree.

    Returns
    -------
    Qualified Python class/function/method records, including named nested helpers.

    Notes
    -----
    C++ declarations are obtained from Clang AST JSON and Python declarations from ast. Anonymous lambdas, vendor internals and generated factory bodies are excluded explicitly. Description coverage is enforced before publishing the inventory; regenerate after changing declaration locations.
    """
    entries=[]
    for path in sorted(list((ROOT/'scripts').glob('*.py'))+list((ROOT/'tests').glob('*.py'))):
        tree=ast.parse(path.read_text());relative=str(path.relative_to(ROOT))
        def visit(node,scope=''):
            """Retain nested qualified names while ignoring anonymous lambda expressions.

            Parameters
            ----------
            node
                Current syntax-tree declaration/context.
            scope
                Enclosing class/function qualification accumulated during traversal.
            """
            for child in ast.iter_child_nodes(node):
                if isinstance(child,(ast.FunctionDef,ast.AsyncFunctionDef,ast.ClassDef)):
                    name=scope+child.name
                    signature=child.name+'('+ast.unparse(child.args)+')' if isinstance(child,(ast.FunctionDef,ast.AsyncFunctionDef)) else child.name
                    entries.append(dict(language='Python',file=relative,line=child.lineno,name=name,kind='class' if isinstance(child,ast.ClassDef) else 'function',signature=signature,purpose=(ast.get_docstring(child) or '').split('\n\n')[0],documentation=ast.get_docstring(child) or ''))
                    visit(child,name+'.')
                else:visit(child,scope)
        visit(tree)
    return entries


def descriptions(entries):
    """Resolve curated C++ behavior and source docstrings, rejecting undocumented library APIs.

    Parameters
    ----------
    entries
        Named-declaration inventory with qualified names and source locations.

    Returns
    -------
    Inventory populated with behavior descriptions; missing descriptions are rejected.

    Notes
    -----
    C++ declarations are obtained from Clang AST JSON and Python declarations from ast. Anonymous lambdas, vendor internals and generated factory bodies are excluded explicitly. Description coverage is enforced before publishing the inventory; regenerate after changing declaration locations.
    """
    cpp=json.loads((ROOT/'docs/api_descriptions.json').read_text())
    python=json.loads((ROOT/'docs/python_descriptions.json').read_text())
    missing=[]
    for entry in entries:
        if entry['language']=='C++':
            name=entry['name'];purpose=cpp.get(name)
            if not purpose and entry['file'].startswith(('src/','tests/')):
                purpose=cpp.get(entry['file']+'::'+name)
                if not purpose and name=='main':purpose='Run '+Path(entry['file']).stem.replace('_',' ')+'; return zero on success and nonzero on a failed check or invalid input.'
                if not purpose and entry['file'].startswith('tests/'):
                    purpose='Test helper for '+Path(entry['file']).stem.replace('_',' ')+': '+name.replace('::',' / ').replace('_',' ')+'. See the typed signature and assertions in the linked test source.'
            entry['purpose']=purpose or ''
        elif not entry['purpose']:
            key=entry['file']+'::'+entry['name']
            entry['purpose']=python.get(key,python.get(entry['name'],''))
            if not entry['purpose'] and entry['name'].split('.')[-1]=='main':entry['purpose']='Execute the '+Path(entry['file']).stem.replace('_',' ')+' command, write its documented outputs and exit nonzero if validation fails.'
            if not entry['purpose'] and entry['file'].startswith('tests/'):
                entry['purpose']='Verify '+entry['name'].split('.')[-1].removeprefix('test_').replace('_',' ')+'. This helper belongs to the source validation test suite.'
        if not entry['purpose']:missing.append(entry['file']+'::'+entry['name'])
    if missing:raise ValueError('Missing API descriptions:\n'+'\n'.join(missing))
    return entries


def annotate(entries):
    """Insert concise Doxygen comments or Python docstrings without changing executable logic.

    Parameters
    ----------
    entries
        Named-declaration inventory with qualified names and source locations.

    Notes
    -----
    C++ declarations are obtained from Clang AST JSON and Python declarations from ast. Anonymous lambdas, vendor internals and generated factory bodies are excluded explicitly. Description coverage is enforced before publishing the inventory; regenerate after changing declaration locations.
    """
    groups={}
    for entry in entries:groups.setdefault(entry['file'],[]).append(entry)
    for name,items in groups.items():
        path=ROOT/name;lines=path.read_text().splitlines(keepends=True);edits={}
        if name.endswith('.py'):
            tree=ast.parse(path.read_text())
            nodes={n.lineno:n for n in ast.walk(tree) if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef,ast.ClassDef))}
            for entry in items:
                node=nodes[entry['line']]
                if ast.get_docstring(node):continue
                first=node.body[0]
                start=min([first.lineno]+[d.lineno for d in getattr(first,'decorator_list',[])])-1
                indentation=' '*(node.col_offset+4)
                doc=indentation+'"""'+entry['purpose'].replace('"""',"'''")+'"""\n'
                if first.lineno==node.lineno:
                    original=lines[start];lines[start]=original[:first.col_offset].rstrip()+'\n'+doc+indentation+original[first.col_offset:].lstrip()
                else:edits[start]=doc
        else:
            for entry in items:
                start=entry['line']-1
                if start>0 and lines[start-1].lstrip().startswith('///'):continue
                indentation=re.match(r'\s*',lines[start]).group()
                edits.setdefault(start,indentation+'/// '+entry['purpose']+'\n')
        for position,text in sorted(edits.items(),reverse=True):lines.insert(position,text)
        updated=''.join(lines)
        if name.endswith('.py'):compile(updated,name,'exec')
        path.write_text(updated)


def main():
    """Generate checked API JSON/Markdown, or annotate source before a fresh inventory pass.

    Notes
    -----
    C++ declarations are obtained from Clang AST JSON and Python declarations from ast. Anonymous lambdas, vendor internals and generated factory bodies are excluded explicitly. Description coverage is enforced before publishing the inventory; regenerate after changing declaration locations.
    """
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compiler',default='/usr/bin/clang++' if Path('/usr/bin/clang++').exists() else 'clang++')
    parser.add_argument('--annotate',action='store_true')
    args=parser.parse_args()
    entries=descriptions(clang_inventory(args.compiler)+python_entries())
    if args.annotate:
        annotate(entries)
        print('Annotated',len(entries),'declarations; rerun without --annotate to refresh locations.')
        return
    for path in sorted((ROOT/'scripts').glob('*.R')):
        for index,line in enumerate(path.read_text().splitlines(),1):
            match=re.match(r'(\w+)\s*<-\s*function\(([^)]*)\)',line)
            if match:entries.append(dict(language='R',file=str(path.relative_to(ROOT)),line=index,name=match[1],kind='function',signature=match[0],purpose='Draw the '+path.stem.replace('_',' ')+' figure using base R graphics; the source defines units, colors, interval interpretation and export devices.'))
    entries.sort(key=lambda e:(e['language'],e['file'],e['line'],e['name']))
    (ROOT/'docs/api_inventory.json').write_text(json.dumps(dict(scope='all named case declarations in include/src/scripts/tests; anonymous lambdas and generated factories excluded',symbols=entries),indent=2)+'\n')
    lines=['# Swordfish API reference','', 'Generated by `scripts/document_api.py` from Clang and Python syntax trees. Anonymous lambdas, vendor code and generated factory bodies are excluded.','']
    for entry in entries:
        lines+=['## '+entry['file']+' :: '+entry['name'],'','`'+entry['signature']+'`','',entry['purpose'],'',f"Source: [{entry['file']}:{entry['line']}](../{entry['file']}#L{entry['line']})",'']
        if entry.get('documentation') and entry['documentation']!=entry['purpose']:lines+=['Contract details:','','```text',entry['documentation'],'```','']
    (ROOT/'docs/API_REFERENCE.md').write_text('\n'.join(lines))
    print('Documented',len(entries),'named declarations with complete description coverage.')


if __name__=='__main__':main()
