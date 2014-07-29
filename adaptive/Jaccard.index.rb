### calculate Jaccard index from exploded files
## already sorted by freq

# column 1 is the NT sequence, column 4 is the freq 


require 'getoptlong'


def main
  minFreq = 0.000001
  n = 500

  optHash = getopt()
  l1 = optHash["--list1"]
  l2 = optHash["--list2"]


  if optHash.key?("--minFreq")
    minFreq = optHash["--minFreq"].to_f
  end
  
  if optHash.key?("--nclone")
    n = optHash["--nclone"].to_i
  end

#  $stderr.puts "#{n}\t#{minFreq}\t#{optHash["--minFreq"]}"

  h1 = readFreq(l1)
  h2 = readFreq(l2)
  e1 = entropy(h1.values)
  e2 = entropy(h2.values)
  c1 = 1 - e1 / Math.log2(h1.size)
  c2 = 1 - e2 / Math.log2(h2.size)
  jsd = jensen_shannon(h1, h2)

  n = [n, h1.size, h2.size].min
  
  a1 = h1.select {|k,v | v >= minFreq}.sort_by {|k,v| v}.reverse[0..n-1]
  a2 = h2.select {|k,v | v >= minFreq}.sort_by {|k,v| v}.reverse[0..n-1]
  
  min1 = a1[-1][1]
  min2 = a2[-1][1]
  
  overlap = findOverlap(a1, a2)
  #  $stderr.puts a1.size
  
  union =  a1.size + a2.size - overlap
  ji = overlap.to_f / union


  puts "#{l1}\t#{l2}\t#{ji.round(5)}\t#{overlap}\t#{union}\t#{min1.round(8)}\t#{min2.round(8)}\t#{e1.round(5)}\t#{e2.round(5)}\t#{c1.round(5)}\t#{c2.round(5)}\t#{jsd.round(5)}"

end

def findOverlap(a1, a2)
  x = a1.to_h.keys
  y = a2.to_h.keys

  return (x & y).size

end

def entropy(f)
  # normalize and only take positives
  f = f.select {|i| i>0}
  s = f.reduce(0) {|sum, x| sum + x}.to_f
  f.map! {|i| i/s}
  return f.reduce(0) {|sum, x| sum + Math.log2(x)*x} * -1
end

def jensen_shannon(a1,a2)
  keys = a1.keys | a2.keys
  hp = entropy(a1.values)
  hq = entropy(a2.values)

  avg = []
  keys.each do |k|
    if !a1.key?(k)
      a1[k] = 0
    end
    if !a2.key?(k)
      a2[k] = 0
    end
    
    avg << 0.5*(a1[k] + a2[k])
  end

  hj = entropy(avg)
  return hj - 0.5*(hp+hq)
end

def readFreq(file)
  h = {}
  t = 0
  File.new(file, "r").each do |line|
    next if line =~ /^nucleotide/
    cols = line.split(/\s+/)
    k , f = cols[0], cols[3].to_f 
    h[k] = f
    t += f
  end
  
  h.each do |k,f|
    f = f.to_f / t
    h[k] = f
  end
  
  return h
end

def getopt
  
  opts = GetoptLong.new(
                        ["--list1", "-a", GetoptLong::REQUIRED_ARGUMENT],
                        ["--list2", "-b", GetoptLong::REQUIRED_ARGUMENT],
                        ["--minFreq", "-f", GetoptLong::OPTIONAL_ARGUMENT],
                        ["--nclone", "-n", GetoptLong::OPTIONAL_ARGUMENT],
                        ["--help", "-h", GetoptLong::NO_ARGUMENT]
                        )
  optHash = {}
  opts.each do |opt, arg|
    optHash[opt] = arg
  end
  
#  $stderr.puts optHash
  
  if optHash.key?("--help") or !optHash.key?("--list1") or !optHash.key?("--list2") 
    $stderr.puts "Usage: ruby __.rb -a sample1 -b sample2 [-f minFreq] [-n n_clones]"
    exit
  end
  return optHash
  
end

main()
