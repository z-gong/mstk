from .ffterm import *
from .forcefield import ForceField

__all__ = ['dff_fuzzy_match']


def dff_fuzzy_match(term, ff):
    '''
    Match a approximate bonded term from a force field by the atom type string.

    It is a implementation of the `parameter transfer` in DFF.

    Parameters
    ----------
    term : BondTerm or AngleTerm or DihedralTerm or ImproperTerm or ChargeIncrementTerm
    ff : ForceField

    Returns
    -------
    matched_term : the same as term
    score : int
    '''
    max_score = 5
    best_score = 0
    best_match = None
    term_candidates = {
        BondTerm           : ff.bond_terms,
        AngleTerm          : ff.angle_terms,
        DihedralTerm       : ff.dihedral_terms,
        ImproperTerm       : ff.improper_terms,
        ChargeIncrementTerm: ff.qinc_terms,
    }
    for term_type, candidates in term_candidates.items():
        if type(term) is term_type:
            for candidate in candidates.values():
                score = _dff_fuzzy_compare(term, candidate, max_score)
                if score > best_score:
                    best_score = score
                    best_match = candidate
            return best_match, best_score

    raise Exception('Invalid type for term')


def _dff_fuzzy_compare(term, candidate, max_score):
    if type(term) is ChargeIncrementTerm and isinstance(candidate, ChargeIncrementTerm):
        score1 = _dff_fuzzy_score(term.type1, candidate.type1, max_score, 3)
        score2 = _dff_fuzzy_score(term.type2, candidate.type2, max_score, 3)
        if score1 * score2 == 0:
            return 0
        return score1 + score2
    if type(term) is BondTerm and isinstance(candidate, BondTerm):
        score1 = _dff_fuzzy_score(term.type1, candidate.type1, max_score, 3)
        score2 = _dff_fuzzy_score(term.type2, candidate.type2, max_score, 3)
        if score1 * score2 == 0:
            return 0
        return score1 + score2
    elif type(term) is AngleTerm and isinstance(candidate, AngleTerm):
        score1 = _dff_fuzzy_score(term.type1, candidate.type1, max_score, 2)
        score2 = _dff_fuzzy_score(term.type2, candidate.type2, max_score, 3)
        score3 = _dff_fuzzy_score(term.type3, candidate.type3, max_score, 2)
        if score1 * score2 * score3 == 0:
            return 0
        return score1 + score2 + score3
    elif type(term) is DihedralTerm and isinstance(candidate, DihedralTerm):
        score1 = _dff_fuzzy_score(term.type1, candidate.type1, max_score, 2)
        score2 = _dff_fuzzy_score(term.type2, candidate.type2, max_score, 3)
        score3 = _dff_fuzzy_score(term.type3, candidate.type3, max_score, 3)
        score4 = _dff_fuzzy_score(term.type4, candidate.type4, max_score, 2)
        if score1 == 0 and (term.type1 == '*' or candidate.type1 == '*'):
            score1 = 1
        if score4 == 0 and (term.type4 == '*' or candidate.type4 == '*'):
            score4 = 1
        if score1 * score2 * score3 * score4 == 0:
            return 0
        return score1 + score2 + score3 + score4
    elif type(term) is ImproperTerm and isinstance(candidate, ImproperTerm):
        score1 = _dff_fuzzy_score(term.type1, candidate.type1, max_score, 3)
        score2 = _dff_fuzzy_score(term.type2, candidate.type2, max_score, 2)
        score3 = _dff_fuzzy_score(term.type3, candidate.type3, max_score, 2)
        score4 = _dff_fuzzy_score(term.type4, candidate.type4, max_score, 2)
        if score2 == 0 and (term.type2 == '*' or candidate.type2 == '*'):
            score2 = 1
        if score3 == 0 and (term.type3 == '*' or candidate.type3 == '*'):
            score3 = 1
        if score4 == 0 and (term.type4 == '*' or candidate.type4 == '*'):
            score4 = 1
        if score1 * score2 * score3 * score4 == 0:
            return 0
        return score1 + score2 + score3 + score4


def _dff_fuzzy_score(name1, name2, max_score, min_length):
    len1, len2 = len(name1), len(name2)
    if len1 < min_length or len2 < min_length:
        return 0
    if name1 == name2:
        return max_score
    elif len1 == len2:
        return _dff_fuzzy_score(name1[:-1], name2[:-1], max_score - 2, min_length)
    elif len1 > len2:
        return _dff_fuzzy_score(name1[:len2], name2, max_score - (len1 - len2), min_length)
    else:
        return _dff_fuzzy_score(name1, name2[:len1], max_score - (len2 - len1), min_length)
