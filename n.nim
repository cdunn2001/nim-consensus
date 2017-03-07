import algorithm
import sequtils
import sets
import strutils

proc get_longest_reads(seqs: seq[string], max_n_read, max_cov_aln: int): seq[string] =
    var longest_n_reads = max_n_read
    if max_cov_aln > 0:
        longest_n_reads = 1
        var seed_len = len(seqs[0])
        var read_cov = 0
        for seq in seqs[1 .. seqs.high]:
            if read_cov div seed_len > max_cov_aln:
                break
            longest_n_reads += 1
            read_cov += len(seq)

        longest_n_reads = min(longest_n_reads, max_n_read)

    #assert longest_n_reads <= len(seqs), "{} <= {}".format(longest_n_reads, len(seqs))
    longest_n_reads = min(longest_n_reads, len(seqs)) # unnec in python
    return seqs[0 .. longest_n_reads-1]
proc get_longest_sorted_reads(seqs: seq[string], max_n_read, max_cov_aln: int): seq[string] =
  # allows us to avoid a redundant sort
  # in get_consensus_trimmed()
  var sorted_seqs: seq[string] = @[]
  sorted_seqs.insert(seqs[1 .. seqs.high], 0)
  #seqs[1 .. seqs.high].sort(proc (x,y: string): int = cmp(y.len, x.len))
  sorted_seqs.sort(proc (x,y: string): int = cmp(y.len, x.len))
  sorted_seqs.insert([seqs[0]], 0)
  echo "sortlen", len(sorted_seqs)
  return get_longest_reads(sorted_seqs, max_n_read, max_cov_aln)
type
  Config = tuple
    max_n_read: int
    min_cov_aln: int
    max_cov_aln: int
iterator get_seq_data(config: Config, min_n_read, min_len_aln: int): auto =
  #let min_idt = 0.70
  #let min_cov = 1
  #let max_n_read 20000
  const max_len = 10000
  var seed_len = 0
  var seed_id: string = ""
  var seqs: seq[string] = @[]
  var read_ids: HashSet[string]
  var read_cov: int
  init(read_ids)
  for line in stdin.lines:
    var l = line.strip.split()
    if len(l) != 2: continue
    var read_id = l[0]
    var sequ = l[1]
    if len(sequ) > max_len:
        sequ = sequ[0 .. max_len-1]
    case read_id
    of "+":
      if len(seqs) >= min_n_read and read_cov div seed_len >= config.min_cov_aln:
          echo "len", len(seqs)
          seqs = get_longest_sorted_reads(seqs, config.max_n_read, config.max_cov_aln)
          yield (seqs, seed_id, config)
      seqs = @[]
      read_ids.init
      seed_id = ""
      read_cov = 0
    of "-":
      break
    of "*":
      seqs = @[]
      read_ids.init
      seed_id = ""
      read_cov = 0
    else:
      if len(sequ) >= min_len_aln:
          if len(seqs) == 0:
              seqs.add(sequ) #the "seed"
              seed_len = len(sequ)
              seed_id = read_id
          if not (read_id in read_ids): #avoidng using the same read twice. seed is used again here by design
              seqs.add(sequ) #the "seed"
              read_ids.incl(read_id)
              read_cov += len(sequ)
proc get_consensus_without_trim(inseqs: seq[string], seed_id: string, config: Config): auto =
    #min_cov, K, max_n_read, min_idt, edge_tolerance, trim_size, min_cov_aln, max_cov_aln = config
    var seqs = inseqs
    if len(seqs) > config.max_n_read:
        seqs = get_longest_sorted_reads(seqs, config.max_n_read, config.max_cov_aln)
    #seqs_ptr = seqs # copy
    #for i in countup(0, len(seqs)-1):
    #  log("i len(seq) %d %d"%(i, len(seq)))
    #consensus_data_ptr = falcon.generate_consensus( seqs_ptr, len(seqs), min_cov, K, min_idt )

    var consensus = "ACGT"
    #consensus = string_at(consensus_data_ptr[0].sequence)[:]
    #eff_cov = consensus_data_ptr[0].eff_cov[:len(consensus)]
    #falcon.free_consensus_data( consensus_data_ptr )
    #del seqs_ptr
    return (consensus, seed_id)
proc main() =
  echo "hi"
  let config: Config = (max_n_read: 500, min_cov_aln: 10, max_cov_aln: 0)
  let min_n_read = 10
  let min_len_aln = 0
  for q in get_seq_data(config, min_n_read, min_len_aln):
    var (seqs, seed_id, config) = q
    echo((len(seqs), seed_id, config))
    var cns: string
    (cns, seed_id) = get_consensus_without_trim(seqs, seed_id, config)
    echo cns, seed_id
when isMainModule:
  main()
