import unittest
#importa i files di aola che ti servono

class TestUtils(unittest.TestCase):

    def test_nome_della_fuznione(self):
        result = nome_della_funzione("parametri di test")
        #assertEqual significa che il risultato deve essere identico a quello specificato da te
        self.assertEqual(result, "valore giusto")
        #puoi usare anche
        #assertNotEqual(result, "valore sbagliato")
        #assertTrue(result)
        #assertFalse(result)
        #assertIn(result, "lista valori giusti")
        #assertNotIn(result, "lista valori sbaliati")
        #assertIsNone(result)

    def test_altra_funzione(self):
        result = nome_altra_funzione("parametri di test")
        #altri test..
if __name__ == '__main__':
    unittest.main()