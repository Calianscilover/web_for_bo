"""Data gates, real lookup replay, paired seeds and online-compatible feedback."""
import csv
import io
import json
from pathlib import Path
import tempfile
import time
import unittest
from unittest.mock import patch
from rdkit import Chem
from rdkit.Chem import Descriptors

import server
from offline_benchmark import dataset as ds
from offline_benchmark.runner import aggregate, run


class OfflineFlow(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.previous = server.OUTPUT
        server.OUTPUT = Path(self.temp.name)
        self.client = server.app.test_client()

    def tearDown(self):
        server.OUTPUT = self.previous
        self.temp.cleanup()

    def example(self):
        result = self.client.post('/api/v1/offline-datasets/example', json={})
        self.assertEqual(result.status_code, 201)
        return result.json['dataset_id']

    def formulate(self, identifier, **kw):
        return self.client.post(f'/api/v1/offline-datasets/{identifier}/validate-formulations', json={
            'validation_mode': 'chemistry', 'input_format': 'experiment_recipe', 'min_salts': 1, 'max_salts': 2,
            'require_all_solvents': True, **kw})

    def ready(self):
        identifier = self.example()
        self.assertEqual(self.formulate(identifier).status_code, 200)
        result = self.client.post(f'/api/v1/offline-datasets/{identifier}/target', json={'column':'CRR[%]'})
        self.assertEqual(result.status_code, 200, result.json)
        self.assertEqual(self.client.post(f'/api/v1/offline-datasets/{identifier}/freeze',json={}).status_code, 200)
        return identifier

    def wait(self, identifier):
        for _ in range(240):
            result = self.client.get(f'/api/v1/offline-benchmarks/{identifier}').json
            if result['status'] not in ('queued','running'):
                return result
            time.sleep(.15)
        self.fail('offline worker timed out')

    def test_gate_sequence_and_real_crr_file(self):
        identifier = self.example()
        root = server.OUTPUT/'offline'/'datasets'/identifier
        self.assertEqual(self.client.post(f'/api/v1/offline-datasets/{identifier}/target',json={'column':'CRR[%]'}).status_code,400)
        self.assertEqual(self.client.post(f'/api/v1/offline-datasets/{identifier}/freeze',json={}).status_code,400)
        self.assertEqual(self.client.post('/api/v1/offline-benchmarks',json={'dataset_id':identifier}).status_code,400)
        result = self.formulate(identifier)
        self.assertEqual(result.status_code,200,result.json)
        self.assertEqual(result.json['N'],41)
        fields, catalog = ds.table(root/'formulation_catalog.csv')
        self.assertNotIn('CRR[%]',fields)
        self.assertNotIn('offline_target',fields)
        self.assertEqual(len([x for x in fields if x.startswith('mass_ratio_')]),8)
        self.assertEqual(self.formulate(identifier,max_salts=1).status_code,400)
        self.assertEqual(self.client.get(f'/api/v1/offline-datasets/{identifier}').json['status'],'uploaded')
        self.assertEqual(self.formulate(identifier).status_code,200)
        result=self.client.post(f'/api/v1/offline-datasets/{identifier}/target',json={'column':'CRR[%]','direction':'max'})
        self.assertEqual(result.json['target_report']['best'],85.44)
        self.assertEqual(self.client.post(f'/api/v1/offline-datasets/{identifier}/target',json={'column':'not_present'}).status_code,400)
        self.assertEqual(self.client.post(f'/api/v1/offline-datasets/{identifier}/freeze',json={}).status_code,400)

    def test_equivalent_candidate_catalog_and_slot_conflict(self):
        identifier=self.example(); self.formulate(identifier)
        root=server.OUTPUT/'offline'/'datasets'/identifier
        fields, rows=ds.table(root/'formulation_catalog.csv')
        stream=io.StringIO(); writer=csv.DictWriter(stream,fieldnames=fields+['CRR[%]']);writer.writeheader()
        writer.writerows({**row,'CRR[%]':42+i} for i,row in enumerate(rows))
        uploaded=self.client.post('/api/v1/offline-datasets',json={'csv':stream.getvalue()}).json['dataset_id']
        result=self.formulate(uploaded,input_format='candidate_catalog')
        self.assertEqual(result.status_code,200,result.json)
        other=server.OUTPUT/'offline'/'datasets'/uploaded
        original=ds.table(root/'formulation_catalog.csv')[1]
        reconstructed=ds.table(other/'formulation_catalog.csv')[1]
        self.assertEqual([x['sample_id'] for x in original],[x['sample_id'] for x in reconstructed])
        for first,second in zip(original,reconstructed):
            for field in fields:
                if field.startswith(('mass_ratio_','mole_ratio_')):
                    self.assertAlmostEqual(float(first[field]),float(second[field]),places=11)
                else:
                    self.assertEqual(first[field],second[field])
        active = next(i for i in range(8) if float(rows[0][f'mass_ratio_{i}']) > 0)
        rows[0][f'compound_{active}']=''
        stream=io.StringIO();writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader();writer.writerows(rows)
        bad=self.client.post('/api/v1/offline-datasets',json={'csv':stream.getvalue()}).json['dataset_id']
        self.assertEqual(self.formulate(bad,input_format='candidate_catalog').status_code,400)

    def test_missing_target_requires_explicit_exclusion_and_freeze_is_immutable(self):
        identifier=self.example()
        root=server.OUTPUT/'offline'/'datasets'/identifier
        fields,rows=ds.table(root/'source.csv');rows[0]['CRR[%]']=''
        ds.write_csv(root/'source.csv',fields,rows)
        self.formulate(identifier)
        url=f'/api/v1/offline-datasets/{identifier}/target'
        self.assertEqual(self.client.post(url,json={'column':'CRR[%]'}).status_code,400)
        result=self.client.post(url,json={'column':'CRR[%]','exclude_missing':True})
        self.assertEqual(result.status_code,200,result.json);self.assertEqual(result.json['N'],40)
        self.client.post(f'/api/v1/offline-datasets/{identifier}/freeze',json={})
        self.assertEqual(self.formulate(identifier).status_code,400)
        self.assertEqual(self.client.post(url,json={'column':'CRR[%]'}).status_code,400)
        (root/'truth.csv').write_text('changed')
        self.assertEqual(self.client.post('/api/v1/offline-benchmarks',json={'dataset_id':identifier}).status_code,400)

    def test_generated_pool_downloaded_recipe_online_offline_roundtrip(self):
        """Exercise the real generate→download own-template→fill→import route."""
        from test_integration import CONFIG
        response=self.client.post('/api/v1/designs',json=CONFIG)
        design=response.json['design_id']
        for _ in range(240):
            status=self.client.get(f'/api/v1/designs/{design}').json
            if status['status'] in ('succeeded','failed'):break
            time.sleep(.05)
        self.assertEqual(status['status'],'succeeded',status)
        own=self.client.get(f'/api/v1/designs/{design}/recipe-template?target=Conductivity')
        self.assertEqual(own.status_code,200)
        examples=list(csv.DictReader(io.StringIO(own.data.decode('utf-8-sig'))))
        downloaded=self.client.get(f'/api/v1/designs/{design}/experiment-template?target=Conductivity')
        reader=csv.DictReader(io.StringIO(downloaded.data.decode('utf-8-sig')))
        fields=reader.fieldnames; rows=list(reader)[:20]
        for i,row in enumerate(rows):row.update(experiment_id='',Conductivity=str(1+i/10))
        stream=io.StringIO();writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader();writer.writerows(examples+rows)
        text=stream.getvalue()
        uploaded=self.client.post('/api/v1/offline-datasets',json={'csv':text}).json
        self.assertEqual(uploaded['inspection']['skipped']['examples'],2)
        identifier=uploaded['dataset_id']
        result=self.formulate(identifier,design_id=design)
        self.assertEqual(result.status_code,200,result.json)
        self.assertEqual(result.json['N'],20)
        self.assertEqual(result.json['property_columns'],['Conductivity'])
        online=self.client.post(f'/api/v1/designs/{design}/experiment-recipes',json={'csv':text,'targets':['Conductivity']})
        self.assertEqual(online.status_code,200,online.json)
        self.assertEqual(online.json['recipes'],20)
        offline_rows=ds.table(server.OUTPUT/'offline'/'datasets'/identifier/'formulation_catalog.csv')[1]
        online_rows=ds.table(server.OUTPUT/'designs'/design/'experiment.csv')[1]
        self.assertEqual({r['sample_id'] for r in offline_rows},{r['sample_id'] for r in online_rows})
        target=self.client.post(f'/api/v1/offline-datasets/{identifier}/target',json={'column':'Conductivity'})
        self.assertEqual(target.json['target_report']['best'],2.9)
        mapping=self.client.get(f'/api/v1/offline-datasets/{identifier}/files/formulation_mapping.csv')
        self.assertEqual(mapping.status_code,200)
        self.assertEqual(len(list(csv.DictReader(io.StringIO(mapping.data.decode('utf-8-sig'))))),20)
        mapping.close()
        self.assertEqual(self.client.get(f'/api/v1/offline-datasets/{identifier}/files/source.csv').status_code,400)

    def test_unielf_split_ions_basis_identity_balance_and_targets(self):
        identifier=self.example();self.formulate(identifier)
        root=server.OUTPUT/'offline'/'datasets'/identifier
        original=ds.table(root/'formulation_catalog.csv')[1][:6]
        rows=[]
        for i,row in enumerate(original):
            species=[];lithium=0
            for j in range(8):
                ratio=float(row[f'mole_ratio_{j}']);name=row[f'compound_{j}'];smiles=row[f'smiles_{j}']
                if ratio==0:continue
                fragments=smiles.split('.')
                if '[Li+]' in fragments:
                    species.append((name+'_1',next(x for x in fragments if x!='[Li+]'),ratio));lithium+=ratio
                else:species.append((name,smiles,ratio))
            species.append(('Li+','[Li+]',lithium));total=sum(x[2] for x in species)
            result={'formulation_id':f'INFER-{i}','CRR[%]':str(40+i) if i<4 else '',
                    'CRR[%]_pred':str(50+i),'CRR[%]_pred_std':'.3','CRR[%]_model_devi':'.5','temperature[C]':'25'}
            for j in range(10):
                n,s,v=species[j] if j<len(species) else ('','','')
                result.update({f'compound_{j}':n,f'smiles_{j}':s,f'ratio_{j}':v/total if v!='' else ''})
            rows.append(result)
        def upload_rows(values):
            stream=io.StringIO();w=csv.DictWriter(stream,fieldnames=list(values[0]));w.writeheader();w.writerows(values)
            return self.client.post('/api/v1/offline-datasets',json={'csv':stream.getvalue()}).json
        uploaded=upload_rows(rows);identifier=uploaded['dataset_id']
        self.assertIn('candidate_catalog',uploaded['formats'])
        self.assertEqual(uploaded['inspection']['property_counts']['CRR[%]'],4)
        self.assertEqual(self.formulate(identifier,input_format='candidate_catalog',catalog_basis='').status_code,400)
        self.assertEqual(self.formulate(identifier,input_format='candidate_catalog',catalog_basis='ratio_mass').status_code,400)
        result=self.formulate(identifier,input_format='candidate_catalog',catalog_basis='unielf_mole')
        self.assertEqual(result.status_code,200,result.json)
        self.assertEqual(result.json['property_columns'],['CRR[%]','CRR[%]_pred'])
        reconstructed=ds.table(server.OUTPUT/'offline'/'datasets'/identifier/'formulation_catalog.csv')[1]
        # No feature column represents Li+ or an anion: all salts are reconstructed.
        self.assertNotIn('Li+',result.json['formulation_report']['components'])
        for original_row,group in zip(original,rows):
            expected={original_row[f'compound_{j}']:float(original_row[f'mass_ratio_{j}']) for j in range(8)}
            found=min(reconstructed,key=lambda r:sum(abs(expected.get(r[f'compound_{j}'],0)-float(r[f'mass_ratio_{j}'])) for j in range(len(result.json['formulation_report']['components']))))
            for j in range(len(result.json['formulation_report']['components'])):
                self.assertAlmostEqual(float(found[f'mass_ratio_{j}']),expected[found[f'compound_{j}']],places=10)
        url=f'/api/v1/offline-datasets/{identifier}/target'
        self.assertEqual(self.client.post(url,json={'column':'CRR[%]'}).status_code,400)
        self.assertEqual(self.client.post(url,json={'column':'CRR[%]_pred'}).status_code,400)
        self.assertEqual(self.client.post(url,json={'column':'CRR[%]_pred_std','source':'external_simulation'}).status_code,400)
        target=self.client.post(url,json={'column':'CRR[%]_pred','source':'external_simulation'})
        self.assertEqual(target.json['N'],6)
        self.assertEqual(target.json['target_report']['best'],55)
        rows[0]['ratio_0']=float(rows[0]['ratio_0'])+.01
        lithium_slot=next(j for j in range(10) if rows[0][f'compound_{j}']=='Li+')
        rows[0][f'ratio_{lithium_slot}']-=.01
        bad=upload_rows(rows)['dataset_id']
        rejected=self.formulate(bad,input_format='candidate_catalog',catalog_basis='unielf_mole')
        self.assertEqual(rejected.status_code,400)
        self.assertIn('不平衡',' '.join(rejected.json['details']))

    def test_all_methods_lookup_paired_initialization_and_budget(self):
        identifier=self.ready()
        methods=['random','qlognei','qlogei','qucb','qpi','qsr','qstd']
        result=self.client.post('/api/v1/offline-benchmarks',json={'dataset_id':identifier,'methods':methods,
            'n_init':8,'total_budget':11,'batch_size':2,'repeats':1,'fit_maxiter':10,'mc_samples':16})
        self.assertEqual(result.status_code,202,result.json)
        bid=result.json['benchmark_id']
        status=self.wait(bid)
        self.assertEqual(status['status'],'completed',status)
        result=self.client.get(f'/api/v1/offline-benchmarks/{bid}/results').json
        self.assertEqual(result['completed'],7)
        self.assertEqual(len({tuple(t['initial_ids']) for t in result['trials']}),1)
        root=server.OUTPUT/'offline'/'benchmarks'/bid
        truth={x['sample_id']:float(x['target_value']) for x in ds.table(server.OUTPUT/'offline'/'datasets'/identifier/'truth.csv')[1]}
        for trial in result['trials']:
            self.assertEqual(len(trial['trace']),11)
            self.assertEqual(len({p['sample_id'] for p in trial['trace']}),11)
            self.assertEqual([p['evaluations'] for p in trial['trace'][-3:]],[10,10,11])
            for point in trial['trace']:
                self.assertEqual(point['value'],truth[point['sample_id']])
            child=root/'methods'/trial['method']/'repeat_000'
            self.assertEqual(len(ds.table(child/'experiment.csv')[1]),8)
            if trial['method']!='random':
                audit=ds.read_json(child/'training_audit_1.json')
                self.assertEqual(audit['n_observed'],8)
                self.assertEqual(set(audit['observed_ids']),set(trial['initial_ids']))
                self.assertEqual(audit['input_dim'],8)
                self.assertTrue(all(x.startswith('mass_ratio_') for x in audit['feature_names']))
                later=ds.read_json(child/'training_audit_2.json')
                self.assertEqual(later['n_observed'],10)
                self.assertEqual(set(later['observed_ids']),{p['sample_id'] for p in trial['trace'][:10]})

    def test_mole_ratio_ions_carbonate_solvents_and_variable_presence(self):
        smiles={'EC':'O=C1OCCO1','DMC':'COC(=O)OC','EMC':'CCOC(=O)OC',
                'MA':'O=C(OC)C','1846':'F[P-](F)(F)(F)(F)F','1954':'[Li+]'}
        compositions=[{'EC':.28,'DMC':.637948718,'1846':.041025641,'1954':.041025641},
                      {'EC':.274647887,'DMC':.536217304,'MA':.108651911,'1846':.040241449,'1954':.040241449},
                      {'EC':.3,'EMC':.6,'1846':.05,'1954':.05},
                      {'EC':.25,'EMC':.55,'MA':.1,'1846':.05,'1954':.05}]
        rows=[]
        for i,composition in enumerate(compositions+[compositions[0]]):
            row={'sample_id':str(i),'temperature[C]':'0','conductivity[mS/cm]':str(3+i),
                 'conductivity_std':'.2','conductivity[mS/cm]_pred':str(4+i),
                 'conductivity[mS/cm]_pred_std':'.5','conductivity[mS/cm]_pred_values':''}
            species=list(composition)
            total_mass=sum(v*Descriptors.MolWt(Chem.MolFromSmiles(smiles[n])) for n,v in composition.items())
            for j in range(10):
                name=species[j] if j<len(species) else ''
                row.update({f'compound_{j}':name,f'smiles_{j}':smiles.get(name,''),
                            f'mole_ratio_{j}':str(composition.get(name,0)),
                            f'mass_ratio_{j}':str(composition[name]*Descriptors.MolWt(Chem.MolFromSmiles(smiles[name]))/total_mass) if name else '0'})
            rows.append(row)
        def upload_rows(data):
            stream=io.StringIO();w=csv.DictWriter(stream,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
            return self.client.post('/api/v1/offline-datasets',json={'csv':stream.getvalue()}).json
        uploaded=upload_rows(rows);identifier=uploaded['dataset_id']
        self.assertEqual(uploaded['inspection']['optional_components'],{'DMC':3,'MA':2,'EMC':2})
        rejected=self.formulate(identifier,input_format='candidate_catalog',catalog_basis='mole')
        self.assertEqual(rejected.status_code,400)
        self.assertIn('取消全部溶剂必含',' '.join(rejected.json['details']))
        validated=self.formulate(identifier,input_format='candidate_catalog',catalog_basis='mole',require_all_solvents=False)
        self.assertEqual(validated.status_code,200,validated.json)
        self.assertEqual(validated.json['N'],4)
        self.assertEqual(validated.json['formulation_report']['replicates_merged'],1)
        self.assertEqual(validated.json['formulation_report']['components'],['LiPF6','EC','DMC','EMC','MA'])
        self.assertEqual(validated.json['property_columns'],['conductivity[mS/cm]','conductivity[mS/cm]_pred'])
        result=self.client.post(f'/api/v1/offline-datasets/{identifier}/target',json={'column':'conductivity[mS/cm]'})
        self.assertEqual(result.status_code,200,result.json)
        self.assertEqual(result.json['target_report']['replicate_measurements'],5)
        self.assertEqual(result.json['N'],4)
        # Source mole fractions describe separate ions. Reconstructed salt fractions
        # must agree with species masses, rather than retaining the ion denominator.
        catalog=ds.table(server.OUTPUT/'offline'/'datasets'/identifier/'formulation_catalog.csv')[1]
        first=min(catalog,key=lambda r:abs(float(r['mass_ratio_1'])-float(rows[0]['mass_ratio_0'])))
        salt_mass=sum(float(rows[0][f'mass_ratio_{j}']) for j in (2,3))
        self.assertAlmostEqual(float(first['mass_ratio_0']),salt_mass,places=10)
        rows[0]['mass_ratio_0']=str(float(rows[0]['mass_ratio_0'])+.01)
        bad=upload_rows(rows)['dataset_id']
        rejected=self.formulate(bad,input_format='candidate_catalog',catalog_basis='mole',require_all_solvents=False)
        self.assertEqual(rejected.status_code,400)
        self.assertIn('质量／摩尔比例不一致',' '.join(rejected.json['details']))

    def test_file_composition_without_registry_aligns_smiles_and_slots(self):
        fields=['sample_id','compound_0','smiles_0','mole_ratio_0','compound_1','smiles_1','mole_ratio_1','y']
        rows=[{'sample_id':'a','compound_0':'new_material','smiles_0':'CC#N','mole_ratio_0':'.2',
               'compound_1':'second_material','smiles_1':'CS(C)=O','mole_ratio_1':'.8','y':'1'},
              {'sample_id':'b','compound_0':'second_alias','smiles_0':'O=S(C)C','mole_ratio_0':'.8',
               'compound_1':'first_alias','smiles_1':'CC#N','mole_ratio_1':'.2','y':'3'},
              {'sample_id':'c','compound_0':'new_material','smiles_0':'CC#N','mole_ratio_0':'.4',
               'compound_1':'second_material','smiles_1':'CS(C)=O','mole_ratio_1':'.6','y':'4'},
              {'sample_id':'d','compound_0':'new_material','smiles_0':'CC#N','mole_ratio_0':'.6',
               'compound_1':'second_material','smiles_1':'CS(C)=O','mole_ratio_1':'.4','y':'6'}]
        def upload_rows(data):
            stream=io.StringIO();writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader();writer.writerows(data)
            return self.client.post('/api/v1/offline-datasets',json={'csv':stream.getvalue()}).json['dataset_id']
        # Even inspection and validation with an empty preset registry must work.
        with patch.object(ds,'LIBRARY',[]):
            identifier=upload_rows(rows)
            result=self.formulate(identifier,validation_mode='composition',input_format='candidate_catalog',catalog_basis='mole')
        self.assertEqual(result.status_code,200,result.json)
        self.assertEqual(result.json['N'],3)
        self.assertEqual(result.json['formulation_report']['replicates_merged'],1)
        self.assertEqual(result.json['formulation_report']['constraints']['source'],'composition_only')
        self.assertEqual(result.json['available_feature_bases'],['mass','mole'])
        target=self.client.post(f'/api/v1/offline-datasets/{identifier}/target',json={'column':'y'})
        self.assertEqual(target.status_code,200,target.json)
        root=server.OUTPUT/'offline'/'datasets'/identifier
        truth=ds.table(root/'truth.csv')[1]
        self.assertEqual(sorted(float(r['target_value']) for r in truth),[2,4,6])
        catalog_fields,_=ds.table(root/'pool_catalog.csv')
        self.assertNotIn('y',catalog_fields)
        rows[1]['compound_1']='new_material';rows[1]['smiles_1']='CCO'
        bad=upload_rows(rows)
        result=self.formulate(bad,validation_mode='composition',input_format='candidate_catalog',catalog_basis='mole')
        self.assertEqual(result.status_code,400)
        self.assertIn('对应不同 SMILES',' '.join(result.json['details']))
        rows[1]['smiles_1']='invalid!'
        bad=upload_rows(rows)
        self.assertEqual(self.formulate(bad,validation_mode='composition',input_format='candidate_catalog',catalog_basis='mole').status_code,400)

    def test_name_only_mass_composition_can_run_bo_and_blocks_mole_features(self):
        text='experiment_id,unknown_A_mass_g,unknown_B_mass_g,y\na,1,4,1\nb,2,3,2\nc,3,2,4\nd,4,1,8\n'
        with patch.object(ds,'LIBRARY',[]):
            identifier=self.client.post('/api/v1/offline-datasets',json={'csv':text}).json['dataset_id']
            result=self.client.post(f'/api/v1/offline-datasets/{identifier}/validate-formulations',json={'input_format':'experiment_recipe'})
        self.assertEqual(result.status_code,200,result.json)
        self.assertEqual(result.json['available_feature_bases'],['mass'])
        root=server.OUTPUT/'offline'/'datasets'/identifier
        fields,rows=ds.table(root/'formulation_catalog.csv')
        self.assertFalse(any(f.startswith('mole_ratio_') for f in fields))
        self.assertTrue(all(not row['smiles_0'] for row in rows))
        enriched=self.client.post('/api/v1/offline-datasets',json={'csv':text}).json['dataset_id']
        info={'components':[{'name':'unknown_A','smiles':'CC#N'},{'name':'unknown_B','smiles':'CS(C)=O'}]}
        augmented=self.formulate(enriched,validation_mode='composition',design_config=info)
        self.assertEqual(augmented.status_code,200,augmented.json)
        self.assertEqual(augmented.json['available_feature_bases'],['mass','mole'])
        self.assertEqual(self.client.post(f'/api/v1/offline-datasets/{identifier}/target',json={'column':'y'}).status_code,200)
        self.assertEqual(self.client.post(f'/api/v1/offline-datasets/{identifier}/freeze',json={}).status_code,200)
        config={'dataset_id':identifier,'n_init':2,'total_budget':4,'repeats':1,'methods':['random','qlogei'],
                'fit_maxiter':10,'mc_samples':16}
        rejected=self.client.post('/api/v1/offline-benchmarks',json={**config,'feature_basis':'mole'})
        self.assertEqual(rejected.status_code,400,rejected.json)
        started=self.client.post('/api/v1/offline-benchmarks',json={**config,'feature_basis':'mass'})
        self.assertEqual(started.status_code,202,started.json)
        status=self.wait(started.json['benchmark_id'])
        self.assertEqual(status['status'],'completed',status)
        results=self.client.get(f'/api/v1/offline-benchmarks/{started.json["benchmark_id"]}/results').json
        self.assertEqual(results['completed'],2)
        self.assertTrue(all(t['first_hit_round'] is not None for t in results['trials']))

    def test_composition_preserves_ions_and_rejects_invalid_ratios(self):
        text='sample_id,compound_0,smiles_0,mole_ratio_0,compound_1,smiles_1,mole_ratio_1,compound_2,smiles_2,mole_ratio_2,y\na,water,O,.8,sodium,[Na+],.1,anion,[Cl-],.1,1\nb,water,O,.6,sodium,[Na+],.2,anion,[Cl-],.2,2\nc,water,O,.4,sodium,[Na+],.3,anion,[Cl-],.3,3\n'
        identifier=self.client.post('/api/v1/offline-datasets',json={'csv':text}).json['dataset_id']
        result=self.formulate(identifier,validation_mode='composition',input_format='candidate_catalog',catalog_basis='mole')
        self.assertEqual(result.status_code,200,result.json)
        self.assertEqual(set(result.json['formulation_report']['components']),{'water','sodium','anion'})
        for bad_text in (text.replace(',O,.8,',',O,-.8,'),text.replace(',O,.8,',',O,.9,'),text.replace(',water,O,.8,',',,,.8,')):
            bad=self.client.post('/api/v1/offline-datasets',json={'csv':bad_text}).json['dataset_id']
            self.assertEqual(self.formulate(bad,validation_mode='composition',input_format='candidate_catalog',catalog_basis='mole').status_code,400)

    def test_resume_reconstructs_pending_round_and_does_not_double_charge(self):
        identifier=self.ready()
        response=self.client.post('/api/v1/offline-benchmarks',json={'dataset_id':identifier,'methods':['random'],
            'n_init':4,'total_budget':9,'batch_size':3,'repeats':1,'fit_maxiter':10,'mc_samples':16})
        bid=response.json['benchmark_id'];self.assertEqual(self.wait(bid)['status'],'completed')
        root=server.OUTPUT/'offline'/'benchmarks'/bid
        original=ds.read_json(root/'results.json')['trials'][0]
        child=root/'methods'/'random'/'repeat_000'
        state=ds.read_json(child/'trial.json');state['status']='cancelled';ds.save_json(child/'trial.json',state)
        fields,rows=ds.table(child/'recommendation_2.csv')
        for row in rows:row['offline_target']=''
        ds.write_csv(child/'recommendation_2.csv',fields,rows)
        state=ds.read_json(root/'status.json');state['status']='cancelled';ds.save_json(root/'status.json',state)
        run(root)
        restored=ds.read_json(root/'results.json')['trials'][0]
        self.assertEqual(original['trace'],restored['trace'])
        self.assertEqual(original['curve'],restored['curve'])
        self.assertEqual(len(restored['trace']),9)

    def test_invalid_budget_and_methods_are_rejected(self):
        identifier=self.ready()
        for extra in ({'n_init':8.5},{'total_budget':42},{'methods':['random','random']},{'methods':['unsupported']}):
            result=self.client.post('/api/v1/offline-benchmarks',json={'dataset_id':identifier,**extra})
            self.assertEqual(result.status_code,400,result.json)

    def test_bound_conditions_and_concentration_are_checked_before_targets(self):
        identifier=self.example()
        result=self.formulate(identifier,design_config={'temperature_C':26})
        self.assertEqual(result.status_code,400)
        result=self.formulate(identifier,design_config={'salt_molality_bounds':[0,.01]})
        self.assertEqual(result.status_code,400,result.json)
        state=self.client.get(f'/api/v1/offline-datasets/{identifier}').json
        self.assertGreater(state['formulation_report']['design_invalid_rows'],0)

    def test_minimization_full_budget_and_two_repeats(self):
        identifier=self.example();self.formulate(identifier)
        result=self.client.post(f'/api/v1/offline-datasets/{identifier}/target',json={'column':'CRR[%]','direction':'min'})
        self.assertEqual(result.json['target_report']['best'],14.93)
        self.client.post(f'/api/v1/offline-datasets/{identifier}/freeze',json={})
        result=self.client.post('/api/v1/offline-benchmarks',json={'dataset_id':identifier,'methods':['random'],
            'n_init':4,'total_budget':41,'batch_size':7,'repeats':2})
        bid=result.json['benchmark_id'];self.assertEqual(self.wait(bid)['status'],'completed')
        trials=self.client.get(f'/api/v1/offline-benchmarks/{bid}/results').json['trials']
        self.assertNotEqual(trials[0]['initial_ids'],trials[1]['initial_ids'])
        for trial in trials:
            self.assertEqual(trial['final_regret'],0)
            first=next(p for p in trial['trace'] if p['value']==14.93)
            self.assertEqual(trial['first_hit_evaluations'],first['evaluations'])
            self.assertEqual(trial['first_hit_round'],first['round'])


if __name__=='__main__':unittest.main()
