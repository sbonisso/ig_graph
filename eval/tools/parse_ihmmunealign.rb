
require 'csv'

class ParseiHMMuneAlign

  def initialize(dir_path, base_str, out_csv="#{dir_path}#{base_str}_preds.csv")
    
    #lst = Dir.glob(dir_path + base_str + "*.txt")
    
    f_lst = dir_path + base_str + "*.txt"
    out_f = dir_path + base_str + "_preds.tab"
    `cat #{f_lst} > #{out_f}`
    #out_csv = dir_path + base_str + "_preds.csv"
    write_csv(out_f, out_csv)
  end

  
  def write_csv(out_f, out_csv)
    self.convert_csv(out_f, out_csv)
  end
  
  
  def self.convert_csv(out_f, out_csv, delim=";")
    id_re = Regexp.new(/IGH(V|D|J)(\d|[a-zA-Z]|\*|-)+/)
    CSV.open(out_csv, "w", {:col_sep => "\t"}) do |csv|
      File.open(out_f) do |f|
        f.each_line do |line|
          #ary = line.split(";")
          ary = line.split(delim)
          row = ary[0..3]
          row[0].delete!("\"")
          
          #puts row.to_s
          
          row[1] = row[1].include?("NA") ? "?,?" : id_re.match(row[1])[0]
          row[3] = row[3].include?("NA") ? "?,?" : id_re.match(row[3])[0]
          if row[2].include?("NO_DGENE_ALIGNMENT") || 
              row[2].include?("NA") then
            row[2] = "?,?"
          elsif /IGHD(.+)\// =~ row[2] then
            a = row[2].split(/(\d+)\/(\d+)/)
            row[2] = [[a[0], a[1], a[3]].join('').delete("\""), 
                      [a[0], a[2], a[3]].join('').delete("\""),].join(",")
          else
            row[2] = id_re.match(row[2])[0]
          end          
          #row[1..3] = row[1..3].map{|s| id_re.match(s)[0]}
          csv << row
        end
      end
    end
  end

  
end


if __FILE__ == $0 then

  require_relative 'parse_ihmmunealign'

  #parse = ParseiHMMuneAlign.new("#{ENV['HOME']}/bin/ihmmune-align/", "smAb_data")
  ParseiHMMuneAlign.convert_csv(ARGV[0], ARGV[1], "\t")
  
end
