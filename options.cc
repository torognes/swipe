/*
    SWIPE
    Smith-Waterman database searches with Inter-sequence Parallel Execution

    Copyright (C) 2008-2014 Torbjorn Rognes, University of Oslo, 
    Oslo University Hospital and Sencel Bioinformatics AS

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact: Torbjorn Rognes <torognes@ifi.uio.no>, 
    Department of Informatics, University of Oslo, 
    PO Box 1080 Blindern, NO-0316 Oslo, Norway
*/

#include "swipe.h"
#include <cassert>
#include <cerrno>  // errno, ERANGE
#include <cinttypes>  // PRId64
#include <cmath>  // std::floor, std::isfinite
#include <cstdint>  // std::int64_t
#include <cstdlib>  // std::strtol, std::strtod
#include <limits>
#include <string>

// the options: parsed and checked by args_init(), shown by args_show()

constexpr int max_threads = 256;

auto args_show(Parameters const & parameters) -> void
{
  if (parameters.view == OutputFormat::plain)
  {
    
    if (cpu_feature_ssse3 == 0)
    {
      fprintf(out, "The performance is reduced because this CPU lacks SSSE3.\n\n");
    }
    
    char const * symtypestring[] = { "Nucleotide", "Amino acid", "Translated query", "Translated database", "Both translated", "Sound" };
    
    //      char * viewtypestring[] = { "plain", 0, 0, 0, 0, 0, 0, "xml",
    //			  "tab-separated", "tab-separated with comments" };
    
    fprintf(out, "Database file:     %s\n", parameters.databasename);
    fprintf(out, "Database title:    %s\n", db_gettitle());
    fprintf(out, "Database time:     %s\n", db_gettime());
    
    if (db_ismasked() != 0)
      {
	fprintf(out, "Database size:     %" PRId64 " residues", db_getsymcount_masked());
	fprintf(out, " in %" PRId64 " sequences\n", db_getseqcount_masked());
      }
      else
      {
	fprintf(out, "Database size:     %" PRId64 " residues", db_getsymcount());
	fprintf(out, " in %" PRId64 " sequences\n", db_getseqcount());
      }

      fprintf(out, "Longest db seq:    %ld residues\n", db_getlongest());

      if (parameters.effdbsize > 0)
      {
	fprintf(out, "Effective db size: %" PRId64 "\n", parameters.effdbsize);
      }

      fprintf(out, "Query file name:   %s\n", parameters.queryname);

      long qlen = 0;
      if ((parameters.symtype == SymbolType::blastn) || (parameters.symtype == SymbolType::blastx) || (parameters.symtype == SymbolType::tblastx))
      {
	qlen = query.nt[0].len;
      }
      else
      {
	qlen = query.aa[0].len;
      }

      fprintf(out, "Query length:      %ld residues\n", qlen);

      query_show();

      if (parameters.symtype == SymbolType::blastn)
      {
	fprintf(out, "Query strands:     ");
	switch (parameters.querystrands)
	{
	case QueryStrands::plus:
	  fprintf(out, "Plus");
	  break;
	case QueryStrands::minus:
	  fprintf(out, "Minus");
	  break;
	case QueryStrands::both:
	  fprintf(out, "Plus and minus");
	  break;
	default:
	  break;
	}
	fprintf(out, "\n");
	fprintf(out, "Score matrix:      %ld/%ld\n", parameters.matchscore, parameters.mismatchscore);
      }
      else
      {
	fprintf(out, "Score matrix:      %s\n", parameters.matrixname);
      }

      fprintf(out, "Gap penalty:       %ld+%ldk\n", parameters.gapopen, parameters.gapextend);
      fprintf(out, "Max expect shown:  %-g\n", parameters.expect);
      fprintf(out, "Min score shown:   %ld\n", parameters.minscore);
      fprintf(out, "Max matches shown: %ld\n", parameters.maxmatches);
      fprintf(out, "Alignments shown:  %ld\n", parameters.alignments);
      fprintf(out, "Show gi's:         %ld\n", parameters.show_gis);
      fprintf(out, "Show taxid's:      %ld\n", parameters.show_taxid);
      fprintf(out, "Threads:           %ld\n", parameters.threads);
      fprintf(out, "Symbol type:       %s\n", symtypestring[static_cast<long>(parameters.symtype)]);
      if ((parameters.symtype == SymbolType::blastx) || (parameters.symtype == SymbolType::tblastx))
      {
	fprintf(out, "Query genetic code:%s (%ld)\n", gencode_names[parameters.query_gencode - 1], parameters.query_gencode);
      }
      if ((parameters.symtype == SymbolType::tblastn) || (parameters.symtype == SymbolType::tblastx))
      {
	fprintf(out, "DB genetic code:   %s (%ld)\n", gencode_names[parameters.db_gencode - 1], parameters.db_gencode);
      }

      // fprintf(out, "View:              %s\n", viewtypestring[view]);
      if (parameters.taxidfilename != nullptr)
      {
	fprintf(out, "Taxid filename:    %s\n", parameters.taxidfilename);
      }
      fprintf(out, "\n");
    }
}
  

// anonymous namespace: limit visibility and usage to this translation unit
namespace {

auto args_usage(char const * const program_name) -> void
{
  /* options unused by BLAST: chkuxHN */
  /* options used by SWIPE:   chkuxHN  */

  fprintf(out, "Usage: %s [OPTIONS]\n", program_name);
  fprintf(out, "  -h, --help                 show help\n");
  fprintf(out, "      --version              show version\n");
  fprintf(out, "  -d, --db=FILE              sequence database base name (required)\n");
  fprintf(out, "  -i, --query=FILE           query sequence filename (stdin)\n");
  fprintf(out, "  -M, --matrix=NAME/FILE     score matrix name or filename (BLOSUM62)\n");
  fprintf(out, "  -q, --penalty=NUM          penalty for nucleotide mismatch (-3)\n");
  fprintf(out, "  -r, --reward=NUM           reward for nucleotide match (1)\n");
  fprintf(out, "  -G, --gapopen=NUM          gap open penalty (11)\n");
  fprintf(out, "  -E, --gapextend=NUM        gap extension penalty (1)\n");
  fprintf(out, "  -v, --num_descriptions=NUM sequence descriptions to show (250)\n");
  fprintf(out, "  -b, --num_alignments=NUM   sequence alignments to show (100)\n");
  fprintf(out, "  -e, --evalue=REAL          maximum expect value of sequences to show (10.0)\n");
  fprintf(out, "  -k, --minevalue=REAL       minimum expect value of sequences to show (0.0)\n");
  fprintf(out, "  -c, --min_score=NUM        minimum score of sequences to show (1)\n");
  fprintf(out, "  -u, --max_score=NUM        maximum score of sequences to show (inf.)\n");
  fprintf(out, "  -a, --num_threads=NUM      number of threads to use [1-%d] (1)\n", max_threads);
  fprintf(out, "  -m, --outfmt=NUM           output format [0,7-9=plain,xml,tsv,tsv+] (0)\n");
  fprintf(out, "  -I, --show_gis             show gi numbers in results (no)\n");
  fprintf(out, "  -p, --symtype=NAME/NUM     symbol type/translation [0-4] (1)\n");
  fprintf(out, "  -S, --strand=NAME/NUM      query strands to search [1-3] (3)\n");
  fprintf(out, "  -Q, --query_gencode=NUM    query genetic code [1-23] (1)\n");
  fprintf(out, "  -D, --db_gencode=NUM       database genetic code [1-23] (1)\n");
  fprintf(out, "  -x, --taxidlist=FILE       taxid list filename (none)\n");
  fprintf(out, "  -N, --dump=NUM             dump database [0-2=no,yes,split headers] (0)\n");
  fprintf(out, "  -H, --show_taxid           show taxid etc in results (no)\n");
  fprintf(out, "  -o, --out=FILE             output file (stdout)\n");
  fprintf(out, "  -z, --dbsize=NUM           set effective database size (0)\n");
}

auto args_version() -> void
{
  char const ref[] = "Reference: T. Rognes (2011) Faster Smith-Waterman database searches\nwith inter-sequence SIMD parallelisation, BMC Bioinformatics, 12:221.";
  fprintf(out, "%s\n\n%s\n", swipe_name_and_version, ref);
}

auto args_help(char const * const program_name) -> void
{
  args_version();
  fprintf(out, "\n");
  
  args_usage(program_name);
}

// strict conversions of option values (KI-9): the whole value must be
// a number, without trailing characters, and within the range of the
// type; otherwise, swipe stops with the error message of the option
auto parse_long(char const * const text, char const * const message) -> long
{
  assert(text != nullptr);
  char * end = nullptr;
  errno = 0;
  auto const value = std::strtol(text, &end, 10);
  if ((end == text) or (*end != '\0') or (errno == ERANGE))
  {
    fatal(message);
  }
  return value;
}

auto parse_double(char const * const text, char const * const message) -> double
{
  assert(text != nullptr);
  char * end = nullptr;
  errno = 0;
  auto const value = std::strtod(text, &end);
  if ((end == text) or (*end != '\0') or (errno == ERANGE) or
      (not std::isfinite(value)))
  {
    fatal(message);
  }
  return value;
}

// the effective database size accepts the real notation of blastall's
// -z (e.g. 7.06e+06, GitHub #9), but must be a non-negative integer
auto parse_dbsize(char const * const text) -> std::int64_t
{
  static char const message[] = "Illegal effective db size specified";
  constexpr auto upper_limit = static_cast<double>(std::numeric_limits<std::int64_t>::max());
  auto const value = parse_double(text, message);
  if ((value < 0.0) or (std::floor(value) < value) or (value >= upper_limit))
  {
    fatal(message);
  }
  return static_cast<std::int64_t>(value);
}

}  // anonymous namespace

auto args_init(int argc, char * const * argv) -> Parameters
{
  Parameters parameters;

  parameters.progname = argv[0];

  opterr = 1;
  char short_options[] = "d:i:M:q:r:G:E:S:v:b:c:u:e:k:a:m:p:x:C:Q:D:F:K:N:o:z:IHh";

  static struct option long_options[] =
  {
    {"db",               required_argument, nullptr, 'd' },
    {"query",            required_argument, nullptr, 'i' },
    {"matrix",           required_argument, nullptr, 'M' },
    {"penalty",          required_argument, nullptr, 'q' },
    {"reward",           required_argument, nullptr, 'r' },
    {"gapopen",          required_argument, nullptr, 'G' },
    {"gapextend",        required_argument, nullptr, 'E' },
    {"strand",           required_argument, nullptr, 'S' },
    {"num_descriptions", required_argument, nullptr, 'v' },
    {"num_alignments",   required_argument, nullptr, 'b' },
    {"min_score",        required_argument, nullptr, 'c' },
    {"max_score",        required_argument, nullptr, 'u' },
    {"evalue",           required_argument, nullptr, 'e' },
    {"minevalue",        required_argument, nullptr, 'k' },
    {"num_threads",      required_argument, nullptr, 'a' },
    {"outfmt",           required_argument, nullptr, 'm' },
    {"symtype",          required_argument, nullptr, 'p' },
    {"taxidlist",        required_argument, nullptr, 'x' },
    {"taxid",            required_argument, nullptr, 'x' },  /* alias (2.1.1 and older) */
    {"comp_based_stats", required_argument, nullptr, 'C' },
    {"query_gencode",    required_argument, nullptr, 'Q' },
    {"db_gencode",       required_argument, nullptr, 'D' },
    {"filter",           required_argument, nullptr, 'F' },
    {"subalignments",    required_argument, nullptr, 'K' },
    {"dump",             required_argument, nullptr, 'N' },
    {"out",              required_argument, nullptr, 'o' },
    {"dbsize",           required_argument, nullptr, 'z' },
    {"show_gis",         no_argument,       nullptr, 'I' },
    {"show_taxid",       no_argument,       nullptr, 'H' },
    {"help",             no_argument,       nullptr, 'h' },
    {"version",          no_argument,       nullptr, 'V' },
    { nullptr, 0, nullptr, 0 },
  };
  
  int option_index = 0;
  int c = 0;

  // gap penalties not given on the command line take the default
  // values of the score matrix or of the symbol type; a penalty of
  // zero is a valid value (KI-6)
  auto gapopen_given = false;
  auto gapextend_given = false;
  
  while (true)
    {
      c = getopt_long(argc, argv, short_options, long_options, &option_index);
      if (c == -1)
      {
	break;
      }

      switch(c)
	{
	case 'a':
	  /* threads */
	  parameters.threads = parse_long(optarg, "Illegal number of threads specified");
	  break;
	  
	case 'b':
	  /* alignments */
	  parameters.alignments = parse_long(optarg, "Illegal number of alignments specified.");
	  break;
	  
	case 'c':
	  /* min score threshold */
	  parameters.minscore = parse_long(optarg, "Illegal minimum score specified.");
	  break;
	  
	case 'C':
	  /* composition-based adjustments */
	  if ((strcasecmp(optarg, "F") != 0) && (strcmp(optarg, "0") != 0))
	  {
	    fatal("Composition-based score adjustments not supported.");
	  }
	  break;

	case 'd':
	  /* database */
	  parameters.databasename = optarg;
	  break;
	  
	case 'D':
	  /* database genetic code */
	  parameters.db_gencode = parse_long(optarg, "Illegal database genetic code specified.");
	  break;
	  
	case 'e':
	  /* evalue */
	  parameters.expect = parse_double(optarg, "Illegal expect value specified.");
	  break;
	  
	case 'E':
	  /* gap extend */
	  parameters.gapextend = parse_long(optarg, "Illegal gap penalties.");
	  gapextend_given = true;
	  break;
	  
	case 'F':
	  /* filter */
	  if ((strlen(optarg) != 0) && (strcasecmp(optarg, "F") != 0))
	  {
	    fatal("Query sequence filtering not supported.");
	  }
	  break;
	  
	case 'G':
	  /* gap open */
	  parameters.gapopen = parse_long(optarg, "Illegal gap penalties.");
	  gapopen_given = true;
	  break;
	  
	case 'h':
	  args_help(parameters.progname);
	  exit(0);
	  break;

	case 'V':
	  /* long option only: -v is --num_descriptions */
	  args_version();
	  exit(0);
	  break;
	  
	case 'H':
	  /* show_taxid */
	  parameters.show_taxid = 1;
	  break;
	  
	case 'i':
	  /* query */
	  parameters.queryname = optarg;
	  break;
	  
	case 'I':
	  /* show_gis */
	  parameters.show_gis = 1;
	  break;
	  
	case 'k':
	  /* min evalue threshold */
	  parameters.minexpect = parse_double(optarg, "Illegal minimum expect value specified.");
	  break;
	  
	case 'K':
	  /* subalignments */
	  parameters.subalignments = parse_long(optarg, "Illegal number of subalignments specified.");
	  break;
	  
	case 'm':
	  /* view */
	  parameters.view = static_cast<OutputFormat>(parse_long(optarg, "Illegal view type."));
	  break;
	  
	case 'M':
	  /* matrix */
	  parameters.matrixname = optarg;
	  break;
	  
	case 'N':
	  /* dump */
	  parameters.dump = parse_long(optarg, "Illegal dump mode.");
	  break;
	  
	case 'o':
	  /* output file */
	  parameters.outfile = optarg;
	  break;
	  
	case 'p':
	  /* symtype */
	  if (strcmp(optarg, "blastn") == 0)
	  {
	    parameters.symtype = SymbolType::blastn;
	  }
	  else if (strcmp(optarg, "blastp") == 0)
	  {
	    parameters.symtype = SymbolType::blastp;
	  }
	  else if (strcmp(optarg, "blastx") == 0)
	  {
	    parameters.symtype = SymbolType::blastx;
	  }
	  else if (strcmp(optarg, "tblastn") == 0)
	  {
	    parameters.symtype = SymbolType::tblastn;
	  }
	  else if (strcmp(optarg, "tblastx") == 0)
	  {
	    parameters.symtype = SymbolType::tblastx;
	  }
	  else if (strcmp(optarg, "sound") == 0)
	  {
	    parameters.symtype = SymbolType::sound;
	  }
	  else
	  {
	    parameters.symtype = static_cast<SymbolType>(parse_long(optarg, "Illegal symbol type."));
	  }
	  break;
	  
	case 'q':
	  /* penalty */
	  parameters.mismatchscore = parse_long(optarg, "Illegal mismatch penalty specified.");
	  break;
	  
	case 'Q':
	  /* query genetic code */
	  parameters.query_gencode = parse_long(optarg, "Illegal query genetic code specified.");
	  break;
	  
	case 'r':
	  /* reward */
	  parameters.matchscore = parse_long(optarg, "Illegal match reward specified.");
	  break;
	  
	case 'S':
	  if (strcmp(optarg, "plus") == 0)
	  {
	    parameters.querystrands = QueryStrands::plus;
	  }
	  else if (strcmp(optarg, "minus") == 0)
	  {
	    parameters.querystrands = QueryStrands::minus;
	  }
	  else if (strcmp(optarg, "both") == 0)
	  {
	    parameters.querystrands = QueryStrands::both;
	  }
	  else
	  {
	    parameters.querystrands = static_cast<QueryStrands>(parse_long(optarg, "Illegal query strands specified."));
	  }
	  break;

	case 'u':
	  /* maxscore */
	  parameters.maxscore = parse_long(optarg, "Illegal maximum score specified.");
	  break;
	  
	case 'v':
	  /* max matches shown */
	  parameters.maxmatches = parse_long(optarg, "Illegal number of descriptions specified.");
	  break;
	  
	case 'x':
	  /* taxid filename */
	  parameters.taxidfilename = optarg;
	  break;
	  
	case 'z':
	  /* effective db size */
	  parameters.effdbsize = parse_dbsize(optarg);
	  break;
	  
	case '?':
	default:
	  args_usage(parameters.progname);
	  exit(1);
	  break;
	}
    }
  
  long gopen_default = 0;
  long gextend_default = 0;

  if (parameters.symtype == SymbolType::blastn)
  {
    if (not gapopen_given)
    {
      parameters.gapopen = 5;
    }
    if (not gapextend_given)
    {
      parameters.gapextend = 2;
    }
  }
  else if (parameters.symtype < SymbolType::sound)
  {
    if (strlen(parameters.matrixname) == 0)
    {
      parameters.matrixname = default_matrixname;
    }

    if (stats_getprefs(parameters.matrixname, & gopen_default, & gextend_default) != 0)
    {
      if (not gapopen_given)
      {
	parameters.gapopen = gopen_default;
      }
      if (not gapextend_given)
      {
	parameters.gapextend = gextend_default;
      }
    }
    else
    {
      // no default for this matrix: a penalty not given is zero
      if ((not gapopen_given) && (not gapextend_given))
      {
	fatal("Unknown score matrix. Gap penalties must be specified (-G and -E).");
      }
    }
  }
  else if (parameters.symtype == SymbolType::sound)
  {
    if (strlen(parameters.matrixname) == 0)
    {
      parameters.matrixname = "IDENTITY_5_1";
    }
    if (not gapopen_given)
    {
      parameters.gapopen = 15;
    }
    if (not gapextend_given)
    {
      parameters.gapextend = 5;
    }
  }

  parameters.gapopenextend = parameters.gapopen + parameters.gapextend;

  if (parameters.effdbsize < 0)
  {
    fatal("Illegal effective db size specified");
  }

  if ((parameters.threads < 1) || (parameters.threads > max_threads))
  {
    fatal("Illegal number of threads specified");
  }

  if (strlen(parameters.databasename) == 0)
  {
    fatal("No database specified.");
  }

  if (!((parameters.view == OutputFormat::plain) || (parameters.view == OutputFormat::xml) || (parameters.view == OutputFormat::tabular) || (parameters.view == OutputFormat::tabular_with_comments) || (parameters.view == OutputFormat::paralign_xml)))
  {
    fatal("Illegal view type.");
  }

  if ((parameters.symtype < SymbolType::blastn) || (parameters.symtype > SymbolType::sound))
  {
    fatal("Illegal symbol type.");
  }

  if ((parameters.gapopen < 0) || (parameters.gapextend < 0) || ((parameters.gapopen + parameters.gapextend) < 1))
  {
    fatal("Illegal gap penalties.");
  }

  if ((parameters.querystrands < QueryStrands::plus) || (parameters.querystrands > QueryStrands::both))
  {
    fatal("Illegal query strands specified.");
  }

  if ((parameters.querystrands == QueryStrands::minus) && ((parameters.symtype == SymbolType::blastp) || (parameters.symtype == SymbolType::tblastn)))
  {
    fatal("Illegal strand specified for protein query.");
  }

  if ((parameters.query_gencode < 1) || (parameters.query_gencode > 23) || (gencode_names[parameters.query_gencode - 1] == nullptr))
  {
    fatal("Illegal query genetic code specified.");
  }

  if ((parameters.db_gencode < 1) || (parameters.db_gencode > 23) || (gencode_names[parameters.db_gencode - 1] == nullptr))
  {
    fatal("Illegal database genetic code specified.");
  }

  if ((parameters.dump < 0) || (parameters.dump > 2))
  {
    fatal("Illegal dump mode.");
  }

  /* ranges of the result limits (KI-7, KI-9) */
  if (parameters.maxmatches < 0)
  {
    fatal("Illegal number of descriptions specified.");
  }

  if (parameters.alignments < 0)
  {
    fatal("Illegal number of alignments specified.");
  }

  /* scores below 1 are not alignments ("Internal error in align
     function.") */
  if (parameters.minscore < 1)
  {
    fatal("Illegal minimum score specified.");
  }

  if (parameters.maxscore < 0)
  {
    fatal("Illegal maximum score specified.");
  }

  if (parameters.expect <= 0.0)
  {
    fatal("Illegal expect value specified.");
  }

  if (parameters.minexpect < 0.0)
  {
    fatal("Illegal minimum expect value specified.");
  }

  /* the output file is opened (and truncated) only once all the
     options are checked (KI-8) */
  if (parameters.outfile != nullptr)
  {
    FILE * f = fopen(parameters.outfile, "w");
    if (f == nullptr)
    {
      fatal("Unable to open output file for writing.");
    }
    out = f;
  }
  
  translate_init(parameters.query_gencode, parameters.db_gencode);

  return parameters;
}
