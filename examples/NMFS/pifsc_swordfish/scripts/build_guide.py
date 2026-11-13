"""Build the swordfish operations and complete named-API PDF guide with ReportLab."""
import argparse
from collections import Counter
from html import escape
import json
import os
from pathlib import Path
import tempfile
import textwrap
from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.platypus import BaseDocTemplate, Frame, PageTemplate, Paragraph, Spacer, PageBreak, Table, TableStyle, KeepTogether, Preformatted, Flowable
from reportlab.platypus.tableofcontents import TableOfContents

ROOT=Path(__file__).resolve().parents[1]
BLUE=colors.HexColor('#173D5E')
TEAL=colors.HexColor('#167D8D')

class GuideTemplate(BaseDocTemplate):
    """Document template with repeatable page furniture, PDF outline and generated contents."""
    def __init__(self,path):
        """Configure US Letter pages, generous margins and the single-column text frame.

        Parameters
        ----------
        path
            Filesystem input/output path in the schema described by the function summary.
        """
        super().__init__(str(path),pagesize=letter,leftMargin=44,rightMargin=44,topMargin=48,bottomMargin=44,title='Quadra Swordfish - Operations and API Guide',author='Quadra swordfish implementation')
        self.addPageTemplates(PageTemplate(id='body',frames=[Frame(self.leftMargin,self.bottomMargin,self.width,self.height,id='body',leftPadding=0,rightPadding=0,topPadding=0,bottomPadding=0)],onPage=self.page_furniture))
    def page_furniture(self,canvas,document):
        """Draw a consistent running title, footer and page number without obscuring content.

        Parameters
        ----------
        canvas
            ReportLab drawing canvas for the current page.
        document
            Active ReportLab document with page/layout state.
        """
        canvas.saveState();canvas.setFillColor(BLUE);canvas.setFont('Helvetica',8)
        canvas.drawString(44,letter[1]-29,'QUADRA / NORTH PACIFIC SWORDFISH')
        canvas.setStrokeColor(colors.HexColor('#D5DFE6'));canvas.line(44,letter[1]-35,letter[0]-44,letter[1]-35)
        canvas.setFillColor(colors.HexColor('#56636B'));canvas.drawString(44,25,'Operations, scientific contracts and named API reference')
        canvas.drawRightString(letter[0]-44,25,str(document.page));canvas.restoreState()
    def afterFlowable(self,flowable):
        """Register section headings in both the contents table and PDF bookmarks.

        Parameters
        ----------
        flowable
            Completed ReportLab flowable inspected for heading navigation.
        """
        if isinstance(flowable,Paragraph) and flowable.style.name in {'Section','Module'}:
            level=0 if flowable.style.name=='Section' else 1;text=flowable.getPlainText()
            key='heading-'+str(self.seq.nextf('heading'))
            self.canv.bookmarkPage(key);self.canv.addOutlineEntry(text,key,level=level,closed=True)
            self.notify('TOCEntry',(level,text,self.page,key))

class Pipeline(Flowable):
    """Vector diagram showing the stage sequence and the publication gate."""
    def __init__(self):
        """Allocate a compact two-row diagram sized to the guide's text frame.
        """
        super().__init__();self.width=524;self.height=103
    def draw(self):
        """Render labeled stages with arrows and a distinct final publish box.
        """
        labels=['Snapshot inputs','Compile + audit','Fit + curvature','Starts + forecasts','Validate + figures','Publish new bundle']
        for index,label in enumerate(labels):
            row=index//3;col=index%3;x=col*178;y=61-row*51
            self.canv.setFillColor(colors.HexColor('#EAF2F6') if index<5 else colors.HexColor('#DDF0E9'))
            self.canv.setStrokeColor(TEAL);self.canv.roundRect(x,y,165,33,5,fill=1,stroke=1)
            self.canv.setFillColor(BLUE);self.canv.setFont('Helvetica-Bold',9);self.canv.drawCentredString(x+82.5,y+13,label)
            if col<2:self.canv.line(x+166,y+16,x+177,y+16)


def build(output,release=None):
    """Create an atomic PDF from the checked API catalog, documented contracts and optional release evidence.

    Parameters
    ----------
    output
        Destination for generated artifacts; the function summary specifies a file or directory.
    release
        Optional accepted bundle providing immutable numerical acceptance evidence.

    Notes
    -----
    The guide is a presentation of the checked API inventory and optional accepted-run evidence. PDF dependencies are separate from the numerical runtime. Rebuilding the PDF does not itself certify model science; render and inspect the resulting pages before delivery.
    """
    inventory=json.loads((ROOT/'docs/api_inventory.json').read_text());entries=inventory['symbols']
    evidence=json.loads((release/'release_manifest.json').read_text()) if release else None
    styles=getSampleStyleSheet()
    styles.add(ParagraphStyle(name='Section',fontName='Helvetica-Bold',fontSize=17,leading=21,textColor=BLUE,spaceBefore=14,spaceAfter=10,keepWithNext=True))
    styles.add(ParagraphStyle(name='Module',fontName='Helvetica-Bold',fontSize=12,leading=16,textColor=TEAL,spaceBefore=12,spaceAfter=7,keepWithNext=True))
    styles.add(ParagraphStyle(name='Symbol',fontName='Helvetica-Bold',fontSize=10,leading=13,textColor=BLUE,spaceBefore=8,spaceAfter=3,keepWithNext=True,wordWrap='CJK'))
    styles.add(ParagraphStyle(name='Text',fontName='Helvetica',fontSize=9,leading=12,spaceAfter=7))
    styles.add(ParagraphStyle(name='SmallNote',fontName='Helvetica',fontSize=7.6,leading=10,textColor=colors.HexColor('#52616B'),spaceAfter=5,wordWrap='CJK'))
    styles.add(ParagraphStyle(name='CodeLine',fontName='Courier',fontSize=7.4,leading=10,spaceAfter=4,keepWithNext=True))
    story=[]
    def paragraph(text,style='Text'):
        """Append escaped prose in a chosen guide style.

        Parameters
        ----------
        text
            Source/document text; output functions preserve intended line breaks.
        style
            Named ReportLab paragraph style.
        """
        story.append(Paragraph(escape(str(text)),styles[style]))
    def heading(text):
        """Begin a top-level guide section and include it in PDF navigation.

        Parameters
        ----------
        text
            Source/document text; output functions preserve intended line breaks.
        """
        paragraph(text,'Section')
    def table(rows,widths):
        """Append a repeating-header table with wrapped cells and restrained grid styling.

        Parameters
        ----------
        rows
            Ordered report records or wrapped document-table rows, depending on the owning function.
        widths
            Document table column widths in PDF points.
        """
        cells=[[Paragraph(escape(str(c)),styles['Text']) for c in r] for r in rows]
        item=Table(cells,colWidths=widths,repeatRows=1,hAlign='LEFT')
        item.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#EAF2F6')),('VALIGN',(0,0),(-1,-1),'TOP'),('LINEBELOW',(0,0),(-1,0),.7,TEAL),('LINEBELOW',(0,1),(-1,-1),.2,colors.HexColor('#D5DFE6')),('LEFTPADDING',(0,0),(-1,-1),7),('RIGHTPADDING',(0,0),(-1,-1),7),('TOPPADDING',(0,0),(-1,-1),6)]))
        story.append(item);story.append(Spacer(1,9))
    def code(text):
        """Wrap long signatures/commands into fixed-width lines that fit the page.

        Parameters
        ----------
        text
            Source/document text; output functions preserve intended line breaks.
        """
        for line in text.splitlines():
            wrapped=textwrap.wrap(line,width=85,break_long_words=True,break_on_hyphens=False) or ['']
            story.append(Preformatted('\n'.join(wrapped),styles['CodeLine']))
    story.append(Spacer(1,55));paragraph('Swordfish model','Section');paragraph('Operations and complete named API guide','Section')
    paragraph('Quadra / PIFSC-ISC North Pacific swordfish / SS3 v3.30.20 comparison')
    story.append(Spacer(1,18));paragraph('Release posture','Module')
    paragraph('Operational hardening and numerical acceptance are documented here. The supported case remains subject to independent scientific review; published bootstrap forecasts have not been reproduced.')
    counts=Counter(e['language'] for e in entries)
    table([['Guide coverage','Evidence'],['Named declarations',str(len(entries))+' ('+', '.join(f'{k}: {v}' for k,v in sorted(counts.items()))+')'],['Scope','Case headers, executables, Python scripts, R plot functions and tests; focused shared-framework appendix'],['Acceptance status',('Passed on '+evidence['environment']['platform']) if evidence else 'Release run not supplied to this document build'],['Source and output identity','SHA-256 fingerprints; clean source copy; new non-overwriting bundle']], [145,379])
    paragraph('Anonymous lambdas, third-party code and generated factory bodies are excluded from named-declaration coverage. Generated factories and their units are described in the guide. No striped marlin code is implied by this swordfish release.','SmallNote')
    story.append(PageBreak());heading('Contents');contents=TableOfContents();contents.levelStyles=[ParagraphStyle(name='Contents0',fontName='Helvetica-Bold',fontSize=9.5,leading=14,leftIndent=0,spaceBefore=5),ParagraphStyle(name='Contents1',fontName='Helvetica',fontSize=8,leading=11,leftIndent=12)];story.append(contents)
    story.append(PageBreak());heading('1. Supported operation and release lifecycle');story.append(Pipeline())
    paragraph('Use scripts/run_release.py for unattended acceptance and result publication. It copies the source and all eight supplied inputs, builds without reusing case binaries or covariance, runs each gate in order, and publishes the complete bundle only after success.')
    code('python3 scripts/run_release.py --output runs/swordfish-release-003')
    table([['Requirement','Contract'],['Runtime','Python 3.9+, C++17 compiler, make and base R on PATH. CXX is a simple executable, not a shell expression.'],['Original inputs','Eight source hashes must match reference/source_manifest.json. Files stay local and are not included in hosted CI.'],['Output','Destination must not exist. New stage logs, frozen workspace and release_manifest.json are published together.'],['Failure','A .failed bundle retains logs and failure.json. No successful destination is published. A lock prevents cooperative writers sharing one destination.']], [110,414])
    paragraph('Atomic publication uses a same-filesystem directory rename. It protects the supported workflow from partial stage results, not hostile filesystem edits, power loss, distributed filesystem behavior, or external writers ignoring the lock. An uncatchable termination can leave a lock requiring operator inspection.')
    paragraph('Low-level make targets are scratch/developer commands and can overwrite build outputs. Their individual files are not completed release bundles. A nonconverged fit now exits before replacing fitted outputs.')
    heading('2. Covariance identity and numerical contracts')
    paragraph('The covariance stamp binds the full native fitted parameter vector, original inputs, configuration, case source/tests/scripts, framework source and bundled numerical dependencies. It hashes both Hessians, the coordinate map, the covariance and curvature diagnostics. Begin/finish identities must match; consumers fail closed if any stamped artifact or source changes.')
    code('make uncertainty\nmake verify-covariance\nmake robustness\nmake projection-uncertainty')
    paragraph('Conservative source hashing includes comments. After any source change, recalculate curvature rather than forcing a timestamp or fabricating provenance. The coordinate reader also verifies coordinate order, bounds, fitted center, matrix symmetry and positive definiteness. No ridge or eigenvalue clipping is applied.')
    paragraph('The likelihood is a joint penalized target with 84 independent optimized coordinates. Recruitment deviations are not Laplace-integrated. Main deviations sum to zero through a dependent final year. Scalar-generic biology and population functions support double and AD values; zero_state explicitly initializes AD numerical zeros.')
    heading('3. Scientific scope, units and conventions')
    table([['Quantity','Convention'],['Population','Two sexes, ages 0-15+, one area; 19 catch fleets plus eight indices. Numbers are thousands of fish.'],['Calendar','1975-2021; four quarters and two subseasons. Age increments occur January 1; settlement is July. Future years are 2022-2031.'],['Weight / biomass','Individual weight is kg; thousands of fish times kg gives tonnes. Growth lengths are cm.'],['Fishing','Apical mortality is an annual rate. Reported annual F pools ages 1-10 under SS3 units 4 and can change during transients even at constant apical rates.'],['Spawning','Legacy v3.30.20 survival timing is deliberately retained for benchmark comparison. A separate calendar-correct helper is explicitly named.'],['Reference values','Original controls and fitted report values are distinct. Fixed male growth L1 is 83.4 in Report.sso versus 83.2 in the control initial value.']], [110,414])
    paragraph('Supported input layouts fail on unimplemented configurations. This code is not a general SS3 interpreter. Parsing retains source locations; comments in source files or attachments provide provenance rather than executable instructions.')
    paragraph('The five forecasts resolve dynamic 20% baseline, mean historical 2008-2010 raw SPR, SPR30, MSY and status quo 2019-2021 targets. Each uses the 2019-2021 fleet/season mixture. Selectivity equality over those years is checked. Counterfactual unfished histories retain deviations and bias adjustment while recomputing BH recruitment under unfished SSB.')
    heading('4. Projection uncertainty and interpretation')
    paragraph('Central differences of the entire forecast mapping at two step sizes propagate the full 84-coordinate covariance as J V J-transpose. Each perturbation rebuilds history, the starting population, fleet mixture and scenario target. This is a recalibrated target rule, not uncertainty conditional on frozen numerical fishing rates.')
    paragraph('Native SEs come from the covariance diagonal. Positive pointwise 95% intervals are estimate times exp(plus/minus 1.959964 times SE/estimate). These local log-scale Wald intervals omit future recruitment process variability, structural uncertainty and bootstrap resampling. They are not simultaneous trajectory bands, realized-outcome prediction intervals, or management-threshold probabilities.')
    paragraph('Published projection differences remain unresolved and can reach about 50%. The source methods and table disagree on historical reference years. The original scenario forecast settings and 100 bootstrap fitted models are required for exact replay; covariance intervals do not resolve those gaps.')
    rows=[['Scenario','2031 SSB, t','Local 95% interval, t']]
    forecast=ROOT/'build/projections/projection_uncertainty.csv'
    if release:forecast=release/'workspace/examples/NMFS/pifsc_swordfish/build/projections/projection_uncertainty.csv'
    if forecast.exists():
        import csv
        for r in csv.DictReader(forecast.open()):
            if r['year']=='2031' and r['quantity']=='spawning_biomass_tonnes':rows.append([r['name'],f"{float(r['estimate']):,.0f}",f"{float(r['lower95']):,.0f} - {float(r['upper95']):,.0f}"])
        table(rows,[180,120,224])
    heading('5. Acceptance evidence and remaining work')
    paragraph('The main make test target now includes projection stationarity, benchmark roots, correlated delta propagation, strict CSV failure paths and runtime publication/provenance tests, in addition to the original SS3 audits and full gradient checks. The release pipeline adds fit, perturbed starts, curvature, uncertainty and figure generation.')
    if evidence:
        table([['Stage','Exit','Seconds']]+[[r['name'],r['returncode'],f"{r['seconds']:.1f}"] for r in evidence['stages']],[290,60,174])
        paragraph('Compiler: '+evidence['environment']['compiler'],'SmallNote')
        paragraph('Python: '+evidence['environment']['python'].splitlines()[0]+'; '+evidence['environment']['R'],'SmallNote')
        paragraph('Published manifest: '+str(release/'release_manifest.json'),'SmallNote')
    paragraph('Ubuntu/macOS portable, framework and sanitizer CI is configured but must be pushed and run before hosted results can be claimed. A locally passed bundle records the platform actually tested. Independent fisheries-model review, retrospectives, broader structural sensitivities, and original published-bootstrap replay remain outstanding scientific work.')
    story.append(PageBreak());heading('6. API reading guide and generated artifacts')
    paragraph('Each API entry below gives the named declaration, source location, typed signature, behavior and public data members where applicable. Units and supported-layout preconditions above apply to all modules. A function returning an aggregate leaves ownership with the caller; explicit mutators such as settle and output writers are identified in their descriptions.')
    table([['Module','Responsibility'],['data / csv / reference_io','Validated schema, indexing, strict numeric/CSV parsing and versioned fixture/native fitted values.'],['dynamics / biology / population','Mortality, settlement/aging, ALKs, selectivity, initialization, hybrid catches and historical recruitment.'],['likelihood / model','Profiled observation predictions, weighted likelihood and bounded optimizer coordinates.'],['benchmarks / projection / scenarios','Equilibrium roots/MSY, deterministic future populations and shared scenario resolution.'],['delta_method / uncertainty executable','Full covariance propagation, interval validation and output correlation matrices.'],['runtime / run_release','Atomic JSON, hashes, covariance identity, isolated builds, gates, logs and publication.'],['prepare / compare / plots / docs','Input factories, SS3 reference extraction, independent output checks, scientific figures and this guide.']], [180,344])
    paragraph('Generated factories imported_swordfish_data, imported_swordfish_biology_data and imported_swordfish_population_data are created by prepare_model/generate_header. Generated report fixtures are comparison targets only; independent dynamics read the imported original catches and controls. Do not edit build headers by hand.')
    for language in ['C++','Python','R']:
        story.append(PageBreak());heading('7. '+language+' declaration reference')
        previous=None
        for entry in [e for e in entries if e['language']==language]:
            if entry['file']!=previous:paragraph(entry['file'],'Module');previous=entry['file']
            paragraph(entry['name'],'Symbol')
            code(entry['signature'])
            paragraph(entry['purpose'])
            if entry.get('fields'):
                paragraph('Members: '+'; '.join(f"{f['name']}: {f['type']}" for f in entry['fields']),'SmallNote')
            paragraph(f"{entry['kind']} | {entry['file']}:{entry['line']}",'SmallNote')
    story.append(PageBreak());heading('8. Shared Quadra changes')
    for name,purpose,path in [
        ('quadra::GradientResult','First-order objective, gradient, coordinate values, norm and report values; no Hessian is constructed.','core/autodiff/model_gradient.hpp'),
        ('quadra::euclidean_norm','Compute the ordinary Euclidean norm used to summarize a coordinate gradient.','core/autodiff/model_gradient.hpp'),
        ('quadra::evaluate_gradient','Create a fresh tape, initialize/evaluate the model, reverse first-order adjoints and extract the coordinate gradient. Empty coordinate vectors fail.','core/autodiff/model_gradient.hpp'),
        ('quadra::solve_random_effects_laplace','Return immediately for an empty random-effect index set. Existing latent-mode behavior remains active when random effects exist.','core/laplace.hpp'),
        ('quadra::laplace_eval_at_u_star_persistent_structured','For no random effects, evaluate the joint first-order target directly and return zero Laplace log-determinant/constant terms, avoiding an unnecessary Hessian.','core/optimizer.hpp')]:
        paragraph(name,'Symbol');paragraph(purpose);paragraph(path,'SmallNote')
    paragraph('Regression coverage includes the new fixed-only optimizer test, scalar-generic model and parameter transform tests, and an existing Laplace optimizer test. The third-order comparison regression also checks all six value-based comparisons after removing duplicate upstream definitions. This is focused regression coverage, not a claim that every Quadra feature has been requalified.')
    heading('9. Regeneration and references')
    code('python3 scripts/document_api.py\npython3 scripts/build_guide.py --release runs/swordfish-release-003')
    paragraph('The inventory uses Clang declaration ASTs for named C++ entities and Python ast for named Python entities. Missing descriptions fail generation. It excludes anonymous lambdas, vendor internals and generated factory bodies; R plot functions are included. Source Doxygen comments/docstrings and docs/API_REFERENCE.md provide editable and searchable companions.')
    paragraph('Primary local references: original SS3 data/control/starter/forecast and Report.sso, ISC 2023 swordfish assessment PDF, source_manifest.json, MODEL_AND_COMPARISON.md, FIT_DIAGNOSTICS.md, PROJECTIONS.md, PROJECTION_UNCERTAINTY.md and RELEASE_READINESS.md. Exact source/output identities and the environment are retained with each accepted bundle.')
    output.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.NamedTemporaryFile(suffix='.pdf',dir=output.parent,delete=False) as stream:temporary=Path(stream.name)
    try:
        GuideTemplate(temporary).multiBuild(story)
        os.replace(temporary,output)
    finally:temporary.unlink(missing_ok=True)
    print('Wrote '+str(output))


def main():
    """Parse optional accepted-release evidence and create the stable guide PDF artifact.

    Notes
    -----
    The guide is a presentation of the checked API inventory and optional accepted-run evidence. PDF dependencies are separate from the numerical runtime. Rebuilding the PDF does not itself certify model science; render and inspect the resulting pages before delivery.
    """
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--release',type=Path)
    parser.add_argument('--output',type=Path,default=ROOT/'output/pdf/swordfish_operations_and_api_guide.pdf')
    args=parser.parse_args();build(args.output,args.release.resolve() if args.release else None)


if __name__=='__main__':main()
