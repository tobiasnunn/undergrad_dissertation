CREATE OR REPLACE TABLE animal_hosts AS
SELECT host, animal_group, domestication_status FROM (
  VALUES
    -- Mammals
    ('Bat', 'Mammal', 'Wild'),
    ('Bos indicus', 'Mammal', 'Domestic'),
    ('Bos taurus', 'Mammal', 'Domestic'),
    ('Bos taurus indicus', 'Mammal', 'Domestic'),
    ('Bovine', 'Mammal', 'Domestic'),
    ('bovine', 'Mammal', 'Domestic'),
    ('Canis', 'Mammal', 'Both'),
    ('Canis lupus', 'Mammal', 'Wild'),
    ('Canis lupus familiaris', 'Mammal', 'Domestic'),
    ('Canine', 'Mammal', 'Domestic'),
    ('canine', 'Mammal', 'Domestic'),
    ('canis', 'Mammal', 'Both'),
    ('Capra hircus', 'Mammal', 'Domestic'),
    ('Cat', 'Mammal', 'Domestic'),
    ('cat', 'Mammal', 'Domestic'),
    ('cattle', 'Mammal', 'Domestic'),
    ('Corynorhinus townsendii', 'Mammal', 'Wild'),
    ('cow', 'Mammal', 'Domestic'),
    ('coyote', 'Mammal', 'Wild'),
    ('dairy cattle', 'Mammal', 'Domestic'),
    ('Delphinidae', 'Mammal', 'Wild'),
    ('Dog', 'Mammal', 'Domestic'),
    ('dog', 'Mammal', 'Domestic'),
    ('dolphin', 'Mammal', 'Wild'),
    ('elephant', 'Mammal', 'Both'),
    ('Eptesicus fuscus', 'Mammal', 'Wild'),
    ('Equus', 'Mammal', 'Both'),
    ('Equus asinus', 'Mammal', 'Domestic'),
    ('Equus caballus', 'Mammal', 'Domestic'),
    ('equidae', 'Mammal', 'Both'),
    ('Felidae', 'Mammal', 'Both'),
    ('Feline', 'Mammal', 'Domestic'),
    ('feline', 'Mammal', 'Domestic'),
    ('Felis catus', 'Mammal', 'Domestic'),
    ('forest musk deer', 'Mammal', 'Wild'),
    ('goat', 'Mammal', 'Domestic'),
    ('Holstein cow', 'Mammal', 'Domestic'),
    ('horse', 'Mammal', 'Domestic'),
    ('Leptonychotes weddellii', 'Mammal', 'Wild'),
    ('Leporidae', 'Mammal', 'Both'),
    ('Macropodidae', 'Mammal', 'Wild'),
    ('Macropus', 'Mammal', 'Wild'),
    ('meat', 'Mammal', 'Domestic'),
    ('mice', 'Mammal', 'Both'),
    ('Mink', 'Mammal', 'Both'),
    ('Mirounga leonina', 'Mammal', 'Wild'),
    ('mouse', 'Mammal', 'Both'),
    ('Mus musculus', 'Mammal', 'Both'),
    ('Mus musculus CF-1', 'Mammal', 'Domestic'),
    ('Myotis evotis', 'Mammal', 'Wild'),
    ('Myotis yumanensis', 'Mammal', 'Wild'),
    ('orangutan', 'Mammal', 'Wild'),
    ('Oryctolagus cuniculus subsp. domesticus', 'Mammal', 'Domestic'),
    ('Oryctolagus cuniculus subsp. domesticus (rabbit)', 'Mammal', 'Domestic'),
    ('Phalangeriformes', 'Mammal', 'Wild'),
    ('Phascolarctos cinereus', 'Mammal', 'Wild'),
    ('pig', 'Mammal', 'Domestic'),
    ('Pinnipedia', 'Mammal', 'Wild'),
    ('Sus scrofa', 'Mammal', 'Both'),
    ('Sus scrofa cristatus', 'Mammal', 'Wild'),
    ('Sus scrofa domesticus', 'Mammal', 'Domestic'),
    ('Tapirus terrestris', 'Mammal', 'Wild'),
    ('Tursiops truncatus (Bottlenose dolphin)', 'Mammal', 'Wild'),
    ('water deer', 'Mammal', 'Wild'),
    
    -- Birds
    ('Bird', 'Bird', 'Both'),
    ('chicken', 'Bird', 'Domestic'),
    ('Columba livia', 'Bird', 'Both'),
    ('Dendrocygna viduata', 'Bird', 'Wild'),
    ('falcons', 'Bird', 'Wild'),
    ('Gallus gallus', 'Bird', 'Both'),
    ('Gallus gallus domesticus', 'Bird', 'Domestic'),
    ('Magellanic penguin', 'Bird', 'Wild'),
    ('Mareca penelope', 'Bird', 'Wild'),
    ('migratory bird', 'Bird', 'Wild'),
    ('Psittaciformes', 'Bird', 'Both'),
    ('swallow', 'Bird', 'Wild'),
    
    -- Fish
    ('Acipenser baerii', 'Fish', 'Both'),
    ('Astronotus ocellatus', 'Fish', 'Both'),
    ('Atlantic salmon', 'Fish', 'Both'),
    ('Bighead carp', 'Fish', 'Both'),
    ('Carassius auratus', 'Fish', 'Domestic'),
    ('Carassius auratus auratus', 'Fish', 'Domestic'),
    ('Cottus cognatus (slimy sculpin)', 'Fish', 'Wild'),
    ('Cunner fish', 'Fish', 'Wild'),
    ('Danio rerio', 'Fish', 'Domestic'),
    ('Dicentrarchus labrax', 'Fish', 'Both'),
    ('Dicologlossa cuneata', 'Fish', 'Wild'),
    ('fish', 'Fish', 'Both'),
    ('Gold fish', 'Fish', 'Domestic'),
    ('Ictalurus punctatus', 'Fish', 'Both'),
    ('Larimichthys crocea', 'Fish', 'Both'),
    ('Maccullochella peelii', 'Fish', 'Wild'),
    ('Maccullochella peelii peelii', 'Fish', 'Wild'),
    ('Myoxocephalus thompsonii (Deepwater Sculpin)', 'Fish', 'Wild'),
    ('Oncorhynchus mykiss', 'Fish', 'Both'),
    ('Oreochromis niloticus', 'Fish', 'Both'),
    ('Oreochromis niloticus (Nile tilapia)', 'Fish', 'Both'),
    ('Paralichthys olivaceus', 'Fish', 'Both'),
    ('Perca fluviatilis', 'Fish', 'Wild'),
    ('Perca fluviatilis L.', 'Fish', 'Wild'),
    ('Plecoglossus altivelis', 'Fish', 'Wild'),
    ('Poecilia sphenops', 'Fish', 'Domestic'),
    ('Procypris rabaudi', 'Fish', 'Wild'),
    ('Pterophyllum altum', 'Fish', 'Domestic'),
    ('rainbow trout', 'Fish', 'Both'),
    ('Salmo coruhensis', 'Fish', 'Wild'),
    ('Salmo salar', 'Fish', 'Both'),
    ('Salvelinus fontinalis', 'Fish', 'Both'),
    ('Scomber scombrus', 'Fish', 'Wild'),
    ('Scophthalmus maximus', 'Fish', 'Both'),
    ('Tetraodon nigroviridis', 'Fish', 'Domestic'),
    ('Tor tambroides', 'Fish', 'Wild'),
    ('yellowfin seabream', 'Fish', 'Wild'),
    
    -- Reptiles
    ('Chrysemys picta', 'Reptile', 'Wild'),
    ('crocodile lizard', 'Reptile', 'Wild'),
    ('Pelodiscus sinensis', 'Reptile', 'Both'),
    ('Physignathus cocincinus', 'Reptile', 'Both'),
    ('Reptilia', 'Reptile', 'Both'),
    ('Sternotherus odoratus', 'Reptile', 'Wild'),
    ('Testudines', 'Reptile', 'Both'),
    
    -- Amphibians
    ('Ambystoma altamirani', 'Amphibian', 'Wild'),
    ('Andrias davidianus', 'Amphibian', 'Wild'),
    ('Boana prasina', 'Amphibian', 'Wild'),
    ('Dugesia japonica', 'Amphibian', 'Wild'),
    ('Hoplobatrachus rugulosus', 'Amphibian', 'Wild'),
    ('Phyllomedusa distincta', 'Amphibian', 'Wild'),
    ('Plethodon cinereus', 'Amphibian', 'Wild'),
    ('Plethodon sp.', 'Amphibian', 'Wild'),
    ('Schmidtea mediterranea', 'Amphibian', 'Wild'),
    ('Scinax fuscovarius', 'Amphibian', 'Wild'),
    ('Scinax fuscovarius (frog)', 'Amphibian', 'Wild'),
    
    -- Insects
    ('Aedes aegypti', 'Insect', 'Wild'),
    ('Anopheles gambiae G3 strain (malaria mosquito)', 'Insect', 'Wild'),
    ('Aphidus ervi', 'Insect', 'Wild'),
    ('Aphodiinae', 'Insect', 'Wild'),
    ('bark bettle', 'Insect', 'Wild'),
    ('Bee', 'Insect', 'Both'),
    ('Bombyx mori (L.)', 'Insect', 'Domestic'),
    ('Cacopsylla pyri', 'Insect', 'Wild'),
    ('Cetonia aurata', 'Insect', 'Wild'),
    ('Curculionidae_Coleoptera', 'Insect', 'Wild'),
    ('Diatraea saccharalis', 'Insect', 'Wild'),
    ('Diptera', 'Insect', 'Wild'),
    ('Dolichoderus sibiricus', 'Insect', 'Wild'),
    ('fly', 'Insect', 'Wild'),
    ('Galeruca laticollis', 'Insect', 'Wild'),
    ('Glossina palpalis palpalis', 'Insect', 'Wild'),
    ('Hemlock woolly adelgid', 'Insect', 'Wild'),
    ('Holotrichia oblita', 'Insect', 'Wild'),
    ('Humpbacked fly', 'Insect', 'Wild'),
    ('Hypothenemus hampei', 'Insect', 'Wild'),
    ('Ips typographus', 'Insect', 'Wild'),
    ('Lagria villosa', 'Insect', 'Wild'),
    ('Lepidoptera', 'Insect', 'Wild'),
    ('Leptinotarsa decemlineata', 'Insect', 'Wild'),
    ('Massospora cicadina', 'Insect', 'Wild'),
    ('Musca domestica', 'Insect', 'Wild'),
    ('Scarabaeidae', 'Insect', 'Wild'),
    ('Scarabaeidae_pupa_Coleoptera', 'Insect', 'Wild'),
    ('silkworm', 'Insect', 'Domestic'),
    ('Tenebrio molitor', 'Insect', 'Both'),
    ('Thaumetopoea processionea', 'Insect', 'Wild'),
    ('Zophobas morio', 'Insect', 'Both'),
    
    -- Crustaceans
    ('Artemia sp.', 'Crustacean', 'Both'),
    ('Cyclops', 'Crustacean', 'Wild'),
    ('Daphnia magna', 'Crustacean', 'Both'),
    ('Ligia exotica', 'Crustacean', 'Wild'),
    ('Macrobrachium rosenbergii', 'Crustacean', 'Both'),
    ('Penaeus vannamei', 'Crustacean', 'Both'),
    ('Porcellio scaber (Woodlouse)', 'Crustacean', 'Wild'),
    ('shrimp', 'Crustacean', 'Both'),
    
    -- Nematodes
    ('Bursaphelenchus xylophilus', 'Nematode', 'Wild'),
    ('Caenorhabditis elegans', 'Nematode', 'Domestic'),
    ('Caenorhabditis elegans MY316', 'Nematode', 'Domestic'),
    ('Caenorhabditis elegans MY379', 'Nematode', 'Domestic'),
    ('Litoditis marina', 'Nematode', 'Wild'),
    ('Meloidogyne sp.', 'Nematode', 'Wild'),
    ('nematode', 'Nematode', 'Both'),
    ('Steinernema sp. 75', 'Nematode', 'Wild'),
    ('Steinernema sp. S86', 'Nematode', 'Wild'),
    
    -- Annelids (worms)
    ('Earthworm (Lumbricus terrestris)', 'Annelid', 'Wild'),
    ('Eisenia andrei', 'Annelid', 'Both'),
    
    -- Myriapods (centipedes/millipedes)
    ('Geophilidae_myriapod', 'Myriapod', 'Wild'),
    ('Lithobius', 'Myriapod', 'Wild'),
    ('Lithobius_myriapod', 'Myriapod', 'Wild'),
    ('Trigoniulus corallinus', 'Myriapod', 'Wild'),
    
    -- Termites/Isoptera
    ('Furculitermes sp.', 'Insect', 'Wild'),
    ('Procubitermes undulans', 'Insect', 'Wild'),
    
    -- Molluscs
    ('marine bivalve', 'Mollusc', 'Wild'),
    ('marine bivalves', 'Mollusc', 'Wild'),
    
    -- Cnidarians (jellyfish, hydra, sea anemones)
    ('Hydra oligactis strain St.-Petersburg', 'Cnidarian', 'Domestic'),
    ('Hydra vulgaris AEP', 'Cnidarian', 'Domestic'),
    ('Nematostella vectensis', 'Cnidarian', 'Domestic'),
    
    -- Sponges
    ('Chondrocladia robertballardi', 'Sponge', 'Wild'),
    ('Cinachyrella kuekenthali', 'Sponge', 'Wild'),
    ('Coscinoderma mathewsi', 'Sponge', 'Wild'),
    ('Dercitus (Stoeba) latex', 'Sponge', 'Wild'),
    ('Lubomirskia baikalensis', 'Sponge', 'Wild'),
    ('Spongilla sp.', 'Sponge', 'Wild'),
    ('Xestospongia sp.', 'Sponge', 'Wild'),
    
    -- Protists
    ('Dictyostelium discoideum', 'Protist', 'Domestic'),
    ('Plagiopyla sp. strain OYSTR', 'Protist', 'Wild'),
    ('ciliate', 'Protist', 'Wild'),
    
    -- Microalgae (also protists/algae)
    ('Asterionella formosa BG1', 'Algae', 'Domestic'),
    ('Chaetoceros calcitrans', 'Algae', 'Both'),
    ('microalgae', 'Algae', 'Both')
    
) AS t(host, animal_group, domestication_status);

CREATE OR REPLACE TABLE raw_pseudomonas_annotations AS
SELECT * FROM read_csv(
    './repos/undergrad_dissertation/02_analysis/post_annotation_analysis/pseudomonas_annotations_kos.tsv.gz',
    delim = '\t',
    header = false,
    columns = {
        'Accession_with_count': 'VARCHAR',
        'KEGG_KO': 'VARCHAR'
    }
);

--clean data and remove rows with no KO value, filter using quality criteria, vertebrate taxa

CREATE OR REPLACE TABLE pseudomonas_annotations AS
SELECT Q1.Accession, Q1.KEGG_ko, ah.animal_group FROM
(SELECT
  LEFT(Accession_with_count, 15) AS Accession, 
  regexp_replace(
    unnest(
      string_split(KEGG_ko, ',')
      ), '^ko:', ''
    ) AS KEGG_ko
FROM raw_pseudomonas_annotations rpa 
WHERE KEGG_KO <> '-') AS Q1
INNER JOIN pseudomonas_metadata pm ON pm.accession = Q1.Accession
INNER JOIN animal_hosts ah ON ah.host = pm.assemblyInfo.biosample.host
WHERE ah.domestication_status IN ('Wild', 'Both') AND
pm.checkmInfo.contamination < 16 AND
pm.checkmInfo.completeness > 79 AND
ah.animal_group IN ('Mammal', 'Bird', 'Amphibian', 'Fish');